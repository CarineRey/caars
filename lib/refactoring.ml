open Core_kernel
open Bistro

let caars_img = Bistro.Shell_dsl.[ docker_image ~account:"carinerey" ~name:"caars_env" ~tag:"master_20200421" () ]

type url = string

type 'a assoc = (string, 'a) List.Assoc.t

type ('a,'b) either =
  | Left of 'a
  | Right of 'b

type single_end_orientation =
  | F
  | R
  | US

type paired_end_orientation =
  | FR
  | RF
  | UP

let ( $ ) a k = List.Assoc.find_exn ~equal:String.equal a k

let assoc keys ~f =
  List.map keys ~f:(fun k -> k, f k)

let assoc_map a ~f =
  List.map a ~f:(fun (k, v) -> k, f k v)

let assoc_filter_map a ~f =
  List.filter_map a ~f:(fun (k, v) -> Option.map (f k v) ~f:(fun v -> k, v))

let assoc_opt keys ~f =
  List.filter_map keys ~f:(fun k ->
      Option.map (f k) ~f:(fun v -> k, v)
    )

let id_concat xs =
  List.filter xs ~f:(String.( <> ) "")
  |> String.concat ~sep:"_"

let descr ?tag d =
  match tag with
  | None -> d
  | Some tag -> sprintf "%s:%s" d tag

class type blast_db = object
  inherit text
  method format : [`blast_db]
end

module OSE_or_PE = struct
  type 'a t =
    | Single_end of {
        reads : 'a ;
        orientation : single_end_orientation ;
      }
    | Paired_end of {
        reads1 : 'a ;
        reads2 : 'a ;
        orientation : paired_end_orientation ;
      }

  let se reads orientation = Single_end { reads ; orientation }
  let pe reads1 reads2 orientation = Paired_end { reads1 ; reads2 ; orientation }
  let map x ~f = match x with
    | Single_end se -> Single_end { se with reads = f se.reads }
    | Paired_end pe -> Paired_end { pe with reads1 = f pe.reads1 ; reads2 = f pe.reads2 }

  let orientation = function
    | Single_end se -> Left se.orientation
    | Paired_end pe -> Right pe.orientation
end

module Sample_source = struct
  type t =
    | Fastq_file of string OSE_or_PE.t
    | Fasta_file of string OSE_or_PE.t

  let orientation = function
    | Fastq_file x
    | Fasta_file x -> OSE_or_PE.orientation x
end

module Misc_workflows = struct
  open Bistro.Shell_dsl

  let concat ?tag xs =
    let descr = descr ?tag "concat" in
    match xs with
    | [] -> Workflow.shell ~descr [ cmd "touch" [ dest ] ]
    | x :: [] -> x
    | fXs ->
      Workflow.shell ~descr [
        cmd "cat" ~stdout:dest [ list dep ~sep:" " fXs ]
      ]

  let fasta_concat = concat
  let fastq_concat = concat

  class type biopython_sequence_index = object
    inherit binary_file
    method format : [`biopython_sequence_index]
  end

  (* Fasta file and its index must be in the same directory due to
     biopython wich retains the relative path between these 2 files.
     A different location is incompatible with the bistro docker usage
     workflow by worflow. To avoid to cp the complete fasta file we
     use a symbolic link. *)
  let build_biopythonindex ?tag (fasta : fasta file) : biopython_sequence_index file =
    Workflow.shell ~version:1 ~descr:(descr ?tag "build_biopythonindex_fasta.py") [
      mkdir_p dest ;
      within_container caars_img (
        and_list [
          cmd "ln" [ string "-s" ; dep fasta ; dest // "seq.fa" ] ;
          cmd "python" [
            file_dump (string Scripts.build_biopythonindex_fasta) ;
            dest // "index" ;
            dest // "seq.fa" ]
        ]
      )
    ]
end

module Cdhitoverlap = struct
  open Bistro.Shell_dsl

  let cdhitoverlap ?tag ?p ?m ?d (fasta : fasta file) : [`cdhit] directory =
    let out = dest // "cluster_rep.fa" in
    Workflow.shell ~version:1 ~descr:(descr ?tag "cdhitlap") [
      mkdir_p dest;
      cmd "cd-hit-lap" ~img:caars_img [
        opt "-i" dep fasta;
        opt "-o" ident out ;
        option ( opt "-p" float ) p;
        option ( opt "-m" float ) m;
        option ( opt "-d" float ) d;
      ]
    ]

  let cluster_rep dir = Workflow.select dir ["cluster_rep.fa"]
  let cluster dir = Workflow.select dir ["cluster_rep.fa.clstr"]
end

module Trinity = struct
  open Bistro.Shell_dsl
  include Trinity

  let single_stranded_or_unstranded = function
    | F -> string "--SS_lib_type F"
    | R -> string "--SS_lib_type R"
    | US -> string ""

  let paired_stranded_or_unstranded = function
    | RF -> string "--SS_lib_type RF"
    | FR -> string "--SS_lib_type FR"
    | UP -> string ""

  let config_trinity_fasta_paired_or_single = function
    | OSE_or_PE.Single_end se ->
      seq ~sep: " " [ string "--single" ; dep se.reads ; single_stranded_or_unstranded se.orientation ]
    | Paired_end pe ->
      seq ~sep: " " [ string "--left" ; dep pe.reads1 ; string "--right" ; dep pe.reads2 ; paired_stranded_or_unstranded pe.orientation ]

  let fasta_read_normalization_get_output ~fasta ~dest=
    let (vars, code) = match fasta with
      | OSE_or_PE.Single_end _ -> (["DEST", dest;
                                "SINGLELINK", string "`readlink single.norm.fa`"],
                               {| mv $SINGLELINK $DEST/"single.norm.fa"|})
      | Paired_end _ -> (["DEST", dest;
                                "LEFTLINK", string "`readlink left.norm.fa`";
                                "RIGHTLINK", string "`readlink right.norm.fa`"],
                               {|echo $LEFTLINK ; mv $LEFTLINK $DEST/"left.norm.fa"; mv $RIGHTLINK $DEST/"right.norm.fa"|})
    in
    Commons.bash_script vars code

  let fasta_read_normalization
    ?(descr = "")
    max_cov
    ~threads
    ?(memory = 1)
    ?(max_memory = 1)
    (fasta : fasta file OSE_or_PE.t)
    : fasta file OSE_or_PE.t =
  let bistro_memory =
    if max_memory > 2 then Int.(min max_memory (memory * 2))
    else 1
  in
  let given_mem =
    if bistro_memory > 2 then Int.(bistro_memory / 2)
    else 1
  in
  (* reserve more memory by bistro than given to normalization tools*)
  let output_dir =
    Workflow.shell ~descr:("fasta_read_normalization" ^ descr) ~version:2 ~np:threads ~mem:(Workflow.int (1024 * bistro_memory)) [
      mkdir_p dest;
      mkdir_p tmp ;
      within_container caars_img (
        and_list [
          cmd "Trinity" [
            string "--no_version_check";
            opt "--max_memory" ident (seq [ string "$((" ; int given_mem ; string " / 1024))G" ]) ;
            opt "--CPU" ident np ;
            string "--just_normalize_reads";
            opt "--normalize_max_read_cov" int max_cov ;
            config_trinity_fasta_paired_or_single fasta ;
            string "--seqType fa" ;
            opt "--output" seq [ ident tmp ; string "/trinity"] ;
          ];
          cd (tmp // "trinity/insilico_read_normalization") ;
          cmd "sh" [ file_dump (fasta_read_normalization_get_output ~fasta ~dest) ];
        ]
      )
    ]
  in
  match fasta with
  | Single_end se -> OSE_or_PE.se (Workflow.select output_dir ["single.norm.fa"]) se.orientation
  | Paired_end pe ->
    OSE_or_PE.pe
      (Workflow.select output_dir ["left.norm.fa"])
      (Workflow.select output_dir ["right.norm.fa"])
      pe.orientation

  let trinity_fasta
      ?tag
      ?full_cleanup
      ?no_normalization
      ~threads
      ?(memory = 1)
      (sample_fasta : fasta file OSE_or_PE.t)
    : fasta file =
    Workflow.shell ~descr:(descr ?tag "Trinity") ~np:threads ~mem:(Workflow.int (1024 * memory)) [
      mkdir_p dest;
      cmd "Trinity" ~img:caars_img [
        string "--no_version_check";
        opt "--max_memory" ident (seq [ string "$((" ; mem ; string " / 1024))G" ]) ;
        opt "--CPU" ident np ;
        option (flag string "--full_cleanup") full_cleanup ;
        option (flag string "--no_normalize_reads") no_normalization ;
        config_trinity_fasta_paired_or_single sample_fasta;
        string "--seqType fa" ;
        opt "--output" seq [ ident dest ; string "/trinity"] ;
      ];
      cmd "sed" [
        string "-re";
        string {|"s/(>[_a-zA-Z0-9]*)( len=[0-9]* path=.*)/\1/"|};
        string "-i";
        seq [ident dest; string "/trinity.Trinity.fasta";];
      ];
    ]
    |> Fn.flip Workflow.select [ "trinity.Trinity.fasta" ]
end

module Rna_sample = struct
  type t = {
    id : string ;
    group_id : string ;
    species : string ;
    reference_species : string list ;
    sample_file : Sample_source.t ;
    run_trinity : bool ;
    run_transdecoder : bool ;
    run_apytram : bool ;
    precomputed_assembly : string option ;
  }
  type 'a assoc = (t * 'a) list
end

module Apytram = struct

  open Bistro.Shell_dsl

  type compressed_read_db = {
    s : Rna_sample.t ;
    concat_fasta : fasta file;
    index_concat_fasta : Misc_workflows.biopython_sequence_index file;
    rep_cluster_fasta : fasta file;
    reformated_cluster : fasta file;
    index_cluster : Misc_workflows.biopython_sequence_index file;
    cluster_rep_blast_db : blast_db file;
  }

  let string_of_db_type = function
    | Left F -> "F"
    | Left R -> "R"
    | Left US -> "single"
    | Right RF -> "RF"
    | Right FR -> "FR"
    | Right UP -> "paired"

  type apytram_output

  let apytram_multi_species
      ?(descr="")
      ?i
      ?evalue
      ?no_best_file
      ?only_best_file
      ?out_by_species
      ?write_even_empty
      ?id
      ?fid
      ?mal
      ?fmal
      ?len
      ?flen
      ?required_coverage
      ?stats
      ?(threads = 1)
      ?(memory = 1)
      ?time_max
      ~query
      ~fam
      (compressed_reads_dbs : compressed_read_db list) : apytram_output directory =

    let memory = match memory with
      | 0 -> 1
      | _ -> memory
    in

    let formated_db_blasts =
      List.map compressed_reads_dbs ~f:(fun db ->
          seq [dep db.cluster_rep_blast_db ; string "/db:"; string db.s.id]
        )
    in
    let db_types =
      List.map compressed_reads_dbs ~f:(fun {s ; _}->
          seq ~sep:":" [string (string_of_db_type (Sample_source.orientation s.Rna_sample.sample_file)); string s.id]
        )
    in
    let formated_fasta =
      List.map compressed_reads_dbs ~f:(fun db ->
          seq [dep db.concat_fasta ; string ":"; string db.s.id]
        )
    in
    let formated_fastaidx =
      List.map compressed_reads_dbs ~f:(fun db ->
          seq [dep db.index_concat_fasta // "index" ; string ":"; string db.s.id]
        )
    in
    let formated_cluster =
      List.map compressed_reads_dbs ~f:(fun db ->
          seq [dep db.reformated_cluster ; string ":"; string db.s.id]
        )
    in
    let formated_clusteridx =
      List.map compressed_reads_dbs ~f:(fun db ->
          seq [dep db.index_cluster // "index" ; string ":"; string db.s.id]
        )
    in


    Workflow.shell  ~version:5 ~descr:("apytram.py" ^ descr) ~np:threads ~mem:(Workflow.int (memory * 1024)) [
      cmd "apytram.py" ~img:caars_img [
        opt "-q" seq [dep query ; string ":"; string fam] ;
        option (opt "-i" int ) i ;
        option (opt "-e" float ) evalue;
        option (opt "-id" float ) id ;
        option (opt "-fid" float ) fid ;
        option (opt "-mal" int ) mal ;
        option (opt "-fmal" float ) fmal ;
        option (opt "-len" float ) len ;
        option (opt "-flen" float ) flen ;
        option (opt "-required_coverage" float ) required_coverage ;
        option (opt "-time_max" int ) time_max ;
        option (flag string "--stats") stats ;
        option (flag string "--no_best_file") no_best_file ;
        option (flag string "--write_even_empty") write_even_empty ;
        option (flag string "--only_best_file") only_best_file ;
        option (flag string "--out_by_species") out_by_species ;
        opt "-memory" ident (seq [ string "$((" ; mem ; string " / 1024))" ]) ;
        opt "-threads" ident np ;
        opt "-d" ident (seq ~sep:"," formated_db_blasts) ;
        opt "-dt" ident (seq ~sep:"," db_types) ;
        opt "-fa" ident (seq ~sep:"," formated_fasta) ;
        opt "-idx" ident (seq ~sep:"," formated_fastaidx) ;
        flag string "--UseIndex" true;
        opt "-clstr" ident (seq ~sep:"," formated_cluster) ;
        opt "-clstridx" ident (seq ~sep:"," formated_clusteridx) ;
        opt "-out" seq [ident dest ; string "/apytram"] ;
        opt "-log" seq [ident dest ; string "/apytram.log"] ;
        opt "-tmp" ident  ( tmp // "apytram_tmp" ) ;
        flag string "--cds" true;
        (*flag string "--keep_tmp" true;
          opt "-tmp" ident  ( dest // "apytram_tmp" ) ;*)
      ]
    ]

  let get_fasta dir ~family_name ~sample_id =
    let filename = sprintf "apytram.%s.%s.fasta" family_name sample_id in
    Workflow.select dir [ filename ]

end

module Sample_sheet = struct

  type ('a,'b) fastX =
    | Fasta of 'a
    | Fastq of 'b

  type merge_criterion =
    | Merge
    | Length
    | Length_complete

  let parse_fastX_path f = match f, Filename.split_extension f with
    | "-", _ -> None
    | x, (_, Some "fa")-> Some (Fasta x)
    | x, (_, Some "fasta") -> Some (Fasta x)
    | x, (_, Some "fq") -> Some (Fastq x)
    | x, (_, Some "fastq") -> Some (Fastq x)
    | x, (_, Some y) -> failwith ({|Syntax error: sample file extension must be in  ["fa", ".fasta","fq","fastq"] (detected extension: |} ^ y  ^ ",  " ^ x ^" )." )
    |  _  -> failwith ({|Syntax error: sample file extension must be in  ["fa", "fasta","fq","fastq"] |})

  let parse_merge_criterion  = function
    | "merge" -> Merge
    | "length" -> Length
    | "length_complete" -> Length_complete
    | _ -> failwith ({| --merge_criterion must be “length“ or “length_complete” or “merge”. “length” means the longest sequence is selected. “length.complete” : means the largest number of complete sites (no gaps). “merge” means that the set of monophyletic sequences is used to build one long “chimera” sequence corresponding to the merging of them.|})

  let parse_orientation id = function
    | "F"  -> Some (Left F)
    | "R"  -> Some (Left R)
    | "US" -> Some (Left US)
    | "RF" -> Some (Right RF)
    | "FR" -> Some (Right FR)
    | "UP" -> Some (Right UP)
    | "-"  -> None
    | _    -> failwith ({|Syntax error in the sample sheet file (sample -> |} ^ id ^ {|: orientation must be in ["F","R","RF","FR","US","UP"] |})

  let str_list_sample_line l =
    let rec str_elements i = function
      | [] -> ""
      | h::t -> "col #" ^ (string_of_int  i) ^ ": " ^ h ^ ";\n" ^ (str_elements (i+1) t)
    in
    "[" ^ (str_elements 1 l) ^ "]"

  let parse_line_fields_of_rna_conf_file = function
    | [ id ; species ; group_id ; ref_species ; path_fastx_single ; path_fastx_left ; path_fastx_right ; orientation ; run_trinity ; path_assembly ; run_apytram] ->
      let run_transdecoder = true in
      let reference_species = List.sort ~compare:String.compare (String.split ~on:',' ref_species) in
      let run_trinity = match run_trinity with
        | "yes" | "Yes" | "y" | "Y" -> true
        | "no" | "No" | "n" | "N" -> false
        | _ -> failwith ({| Syntax error inthe sample sheet file (sample -> |} ^ id ^ {|): run_trinity must be "yes" or "no" |})
      in
      let precomputed_assembly = match path_assembly, run_trinity with
        | "-", _ -> None
        | path, true -> Some path
        | _, false ->
          failwithf {|Incorrect sample sheet file (sample -> %s): you gave a trinity assembly path but run_trinity is false. It is incompatible.|} id ()
      in
      let run_apytram = match run_apytram with
        | "yes" | "Yes" | "y" | "Y" -> true
        | "no" | "No" | "n" | "N" -> false
        | _ -> failwith ({|Syntax error in the sample sheet file (sample -> |} ^ id ^ {|: run_apytram must be "yes" or "no" |})
      in
      let sample_file= match parse_fastX_path path_fastx_single,
                             parse_fastX_path path_fastx_left,
                             parse_fastX_path path_fastx_right,
                             parse_orientation id orientation,
                             run_apytram,
                             run_trinity with
      | None, Some (Fastq _), Some (Fastq _), Some (Right o), _, _ ->
        Sample_source.Fastq_file (OSE_or_PE.pe path_fastx_left path_fastx_right o)

      | None, Some (Fasta _), Some (Fasta _), Some (Right o), _, _ ->
        Sample_source.Fasta_file (OSE_or_PE.pe path_fastx_left path_fastx_right o)

      | Some (Fastq _), None, None, Some (Left o), _, _ ->
        Sample_source.Fastq_file (OSE_or_PE.se path_fastx_single o)

      | Some (Fasta _), None, None, Some (Left o), _, _ ->
        Sample_source.Fasta_file (OSE_or_PE.se path_fastx_single o)

      | None, None, None, None, false, true ->
        Sample_source.Fastq_file (OSE_or_PE.se "-" US)

      | None, Some (Fasta _ ), Some (Fastq _ ), _, _, _
      | None, Some (Fastq _ ), Some (Fasta _ ), _, _, _ ->
        failwithf {|Incorrect sample sheet file (sample -> %s): you provided a fasta file and a fastq file|} id ()

      | None, None, None, _, true, _ ->
        failwithf {|Incorrect sample sheet file (sample -> %s): you didn't give any RNA-seq data, but you asked to run apytram : it is impossible, apytram needs raw RNA-seq data.|}  id ()

      | Some _, None, None, Some (Right _), _, _ ->
        failwithf {|Incorrect sample sheet file (sample -> %s): Incompatible choice. Path for a single-end data but an orientation for a paired-end data.|} id ()

      | None, Some  _, Some _, Some (Left _), _, _ ->
        failwithf {|Incorrect sample sheet file (sample -> %s): Incompatible choice. Paths for paired-end data but an orientation for a single-end data.|} id ()

      | _, _, _, None, _, _ -> failwith ({|Incorrect sample sheet file (sample -> |} ^ id ^ {|): No given orientation.|})

      | _ -> failwithf {|Incorrect sample sheet file (sample -> %s): Incompatible choices. path_fastq_single must be "-" if data are "paired-end" and path_fastq_left and path_fastq_right must be "-" if data are "single-end".|} id ()
      in
      Some { Rna_sample.id ;
             species ;
             group_id ;
             reference_species ;
             sample_file ;
             run_trinity ;
             run_transdecoder ;
             run_apytram ;
             precomputed_assembly ;
           }
    | l ->
      failwithf
        "Syntax error in sample sheet file. There aren't 11 tab delimited columns. Contents must be:\n\n%s\nbut your contents is:\n\n%s"
        (str_list_sample_line ["id->IDx";"species->my_species";"apytram_group->group1";"ref_species->ref_species";"path_fastq_single->single.fastq/-";"path_fastq_left->-/paired_1.fastq";"path_fastq_right->-/paired_2.fastq";"orientation [F,R,FR,RF,UP,US]";"run_trinity-> y/n";"path_assembly->assembly.fa";"run_apytram-> y/n"])
        (str_list_sample_line l) ()

  let read path =
    let samples =
      In_channel.read_lines path
      |> List.tl_exn (* remove the first line*)
      |> List.map ~f:(String.split ~on:'\t')
      |> List.filter_map ~f:parse_line_fields_of_rna_conf_file
    in
    let id_list = List.map samples ~f:(fun s -> s.id) in
    if List.contains_dup id_list ~compare:String.compare then
      failwith {|There are duplicate ids in the first colum of the sample sheet.|}
    else samples

end

module Family = struct
  type t = {
    name : string ;
    id : int;
  }
end

module Dataset = struct
  type merge_criterion =
    | Merge
    | Length
    | Length_complete

  type t = {
    samples : Rna_sample.t list ;
    sample_sheet : string ;
    alignments_dir : string ;
    seq2sp_dir : string ;
    species_tree_file : string ;
    reference_species : string list ;
    all_families : Family.t list ;
    used_families : Family.t list ;
  }

  let samples_with_trinity_run dataset =
    List.filter dataset.samples ~f:(fun s -> s.run_trinity)

  let families_of_alignments_dir alignments_dir =
    Sys.readdir alignments_dir
    |> Array.filter ~f:(fun f -> Filename.check_suffix f ".fa")
    |> Array.map ~f:(fun f -> fst (String.lsplit2_exn f ~on:'.')) (* Il y a obligatoirement un point dans le nom du fichier fasta *)
    |> Array.to_list

  let parse_family_subset path =
    In_channel.read_lines path
    |> List.map ~f:String.strip
    |> List.filter ~f:(String.( <> ) "")
    |> String.Set.of_list

  let make ?family_subset_file ~sample_sheet ~species_tree_file ~alignments_dir ~seq2sp_dir () =
    let samples = Sample_sheet.read sample_sheet in
    let reference_species =
      List.filter_map samples ~f:(fun s ->
          if s.run_apytram || s.run_trinity then Some s.reference_species
          else None
        )
      |> List.concat
    in
    let all_families_noid = families_of_alignments_dir alignments_dir in
    let used_families_noid =
      match family_subset_file with
      | None -> all_families_noid
      | Some path ->
        let subset = parse_family_subset path in
        let all = String.Set.of_list all_families_noid in
        let diff = String.Set.diff subset all in
        if String.Set.is_empty diff then
          String.Set.to_list subset
        else
          let fams = String.Set.to_list diff |> String.concat ~sep:", " in
          failwithf "Families %s don't exist, see in %s" fams path ()
    in
    let family_ids =
      List.foldi all_families_noid ~init:String.Map.empty ~f:(fun i acc name ->
          String.Map.add_exn acc ~key:name ~data:Family.{ name ; id = i + 1 }
        )
    in
    let used_families = List.map used_families_noid ~f:(fun u_f ->
        String.Map.find_exn family_ids u_f
      )
    in
    if Filename.is_relative species_tree_file then
      failwith {|caars needs the absolute path of the species tree.|}
    else if List.is_empty all_families_noid then
      failwith ({|No files with .fa extention in |} ^ alignments_dir)
    else
      {
        samples ;
        sample_sheet ;
        alignments_dir : string ;
        seq2sp_dir : string ;
        species_tree_file ;
        reference_species : string list ;
        used_families ;
        all_families = String.Map.data family_ids ;
      }

  let apytram_reference_species dataset =
    List.filter_map dataset.samples ~f:(fun s ->
        if s.run_apytram then Some s.reference_species
        else None
      )
    |> List.dedup_and_sort ~compare:(List.compare String.compare)

  let apytram_groups dataset =
    List.filter_map dataset.samples ~f:(fun s ->
        if s.run_apytram then Some s.group_id
        else None
      )
    |> List.dedup_and_sort ~compare:String.compare

  let apytram_samples dataset =
    List.filter dataset.samples ~f:(fun s -> s.run_apytram)
end

module Configuration_directory = struct
  open Bistro.Shell_dsl
  type t = [`caars_configuration] directory

  let make ~memory (config : Dataset.t) : t =
    let families_out = dest // "DetectedFamilies.txt" in
    let script = Bistro.Template_dsl.(
        [
          [seq ~sep:"\t" [string "Detected_families"; string "Fam_ID"]];
          List.map config.all_families ~f:(fun fam -> seq ~sep:"\t" [string fam.name; int fam.id ])
        ]
        |> List.concat
        |> seq ~sep:"\n"
      )
    in
    let ali_cmd_list = List.map config.all_families ~f:(fun fam ->
        let ali = Workflow.input (config.alignments_dir ^ "/" ^ fam.name ^ ".fa") in
        cmd "echo" [dep ali]
      )
    in
    let sp2seq_files = Sys.readdir config.seq2sp_dir
                       |> Array.filter ~f:(fun f ->
                           if Filename.check_suffix f ".tsv" then
                             true
                           else
                             false)
                       |> Array.to_list
    in
    let sp2seq_cmd_list = List.map sp2seq_files ~f:(fun sp2seq_file ->
        let sp2seq = Workflow.input (config.seq2sp_dir ^ "/" ^ sp2seq_file) in
        cmd "echo" [dep sp2seq]
      )
    in
    Workflow.shell ~np:1 ~descr:"parse_input" ~version:17 ~mem:(Workflow.int (memory * 1024)) (List.concat [
        [mkdir_p dest;
         cmd "python" ~img:caars_img [
           file_dump (string Scripts.parse_input) ;
           dep (Workflow.input config.sample_sheet) ;
           dep (Workflow.input config.species_tree_file) ;
           dep (Workflow.input config.alignments_dir) ;
           dep (Workflow.input config.seq2sp_dir) ;
           ident dest ;
         ];
         cmd "cp" [ file_dump script; families_out];
        ];
        ali_cmd_list;
        sp2seq_cmd_list;
      ])

  let ref_transcriptome dir species =
    Workflow.select dir ["R_Sp_transcriptomes" ;  species ^ "_transcriptome.fa" ]

  let ref_seq_fam_links dir species =
    Workflow.select dir ["R_Sp_Seq_Fam_links";  species ^ "_Fam_Seq.tsv"  ]

  let ref_fams dir species family =
    Workflow.select dir ["R_Sp_Gene_Families"; species ^ "." ^ family ^ ".fa"]

  let ali_species2seq_links dir family =
    Workflow.select dir ["Alignments_Species2Sequences" ; "alignments." ^  family ^ ".sp2seq.txt" ]

  let ref_blast_dbs_of_configuration_dir (config : Dataset.t) config_dir =
    List.map config.reference_species ~f:(fun ref_species ->
        let fasta = ref_transcriptome config_dir ref_species in
        let parse_seqids = true in
        let dbtype = "nucl" in
        (ref_species, BlastPlus.makeblastdb ~parse_seqids ~dbtype  ("DB_" ^ ref_species) fasta)
      )
end

module Pipeline = struct
  type t = {
    fasta_reads : fasta file OSE_or_PE.t Rna_sample.assoc ;
    normalized_fasta_reads : fasta file OSE_or_PE.t Rna_sample.assoc ;
    trinity_assemblies : fasta file Rna_sample.assoc ;
    trinity_orfs : fasta file Rna_sample.assoc ;
    trinity_assemblies_stats : text file Rna_sample.assoc ;
    trinity_orfs_stats : text file Rna_sample.assoc ;
    trinity_annotated_fams : [`seq_dispatcher] directory Rna_sample.assoc ;
    ref_blast_dbs : blast_db file assoc ;
    reads_blast_dbs : Apytram.compressed_read_db Rna_sample.assoc ;
    apytram_annotated_ref_fams_by_fam_by_groups : (Family.t * (Rna_sample.t * Family.t * fasta file) list) list ;
  }

  let rna_sample_needs_rna (s : Rna_sample.t) =
    match s.run_apytram, s.run_trinity, s.precomputed_assembly with
    | true, _, _          -> true
    | false, true, Some _ -> false
    | false, true, None   -> true
    | false, false, _     -> false

  let fasta_reads (config : Dataset.t) =
    assoc_opt config.samples ~f:(fun s ->
        if rna_sample_needs_rna s then
          let open OSE_or_PE in
          let fq2fa ?(tag = "") x =
            let descr = id_concat [s.id ; s.species ; tag] in
            Trinity.fastq2fasta ~descr x
          in
          match s.sample_file with
          | Fasta_file fa -> Some (OSE_or_PE.map ~f:Workflow.input fa)
          | Fastq_file (Single_end se) ->
            Some (OSE_or_PE.se (fq2fa (Workflow.input se.reads)) se.orientation)
          | Fastq_file (Paired_end pe) ->
            Some (OSE_or_PE.pe
                    (fq2fa ~tag:"left" (Workflow.input pe.reads1))
                    (fq2fa ~tag:"right" (Workflow.input pe.reads2))
                    pe.orientation)
        else None
      )

  let normalize_fasta_reads fasta_reads memory max_memory threads =
     assoc_map fasta_reads ~f:(fun (s : Rna_sample.t) fa ->
        let max_cov = 20 in
        let descr = id_concat [s.id ; s.species] in
        Trinity.fasta_read_normalization ~descr max_cov ~threads ~memory ~max_memory fa
      )

  let trinity_assemblies_of_norm_fasta normalized_fasta_reads ~memory ~nthreads =
    assoc_filter_map normalized_fasta_reads ~f:(fun (s : Rna_sample.t) normalized_fasta_reads ->
        match s.run_trinity, s.precomputed_assembly with
        | true, None ->
          let tag = id_concat [s.id ; s.species] in
          Trinity.trinity_fasta ~tag ~no_normalization:true ~full_cleanup:true ~memory ~threads:nthreads normalized_fasta_reads
          |> Option.some
        | _, Some assembly_path -> Some (Workflow.input assembly_path)
        | (_, _)   -> None
      )

  let transdecoder_orfs_of_trinity_assemblies trinity_assemblies ~memory ~nthreads =
    assoc_map trinity_assemblies ~f:(fun (s : Rna_sample.t) trinity_assembly ->
        match s.run_transdecoder, s.precomputed_assembly with
        | true, None ->
          let pep_min_length = 50 in
          let retain_long_orfs = 150 in
          let descr = id_concat ["Assembly" ; s.id ; s.species] in
          Transdecoder.transdecoder
            ~descr ~retain_long_orfs ~pep_min_length ~only_best_orf:false
            ~memory ~threads:nthreads
            trinity_assembly
        | (false, _ )
        | (true, Some _) -> trinity_assembly
      )

  let assemblies_stats_of_assemblies assemblies =
    assoc_filter_map assemblies ~f:(fun (s : Rna_sample.t) assembly ->
        match s.precomputed_assembly with
        | Some _ -> None
        | None ->
          Some (Trinity.assembly_stats ~descr:(s.id ^ "_" ^ s.species) assembly)
      )

  let ref_blast_dbs (config : Dataset.t) config_dir =
    assoc config.reference_species ~f:(fun ref_species ->
        let fasta = Configuration_directory.ref_transcriptome config_dir ref_species in
        let parse_seqids = true in
        let dbtype = "nucl" in
        BlastPlus.makeblastdb ~parse_seqids ~dbtype  ("DB_" ^ ref_species) fasta
      )

  let reformat_cdhit_cluster ?tag cluster : fasta file =
    let open Bistro.Shell_dsl in
    Workflow.shell ~version:1 ~descr:(descr ?tag "reformat_cdhit_cluster2fasta.py") [
      cmd "python" ~img:caars_img [
        file_dump (string Scripts.reformat_cdhit_cluster2fasta);
        dep cluster ;
        ident dest]
    ]

  let blast_dbs_of_norm_fasta norm_fasta =
    assoc_filter_map norm_fasta ~f:(fun (s : Rna_sample.t) norm_fasta ->
        if s.run_apytram then
          let concat_fasta = match norm_fasta with
            | OSE_or_PE.Single_end se -> se.reads
            | Paired_end pe ->
              Misc_workflows.fasta_concat ~tag:(id_concat [s.id ; "fasta_lr"]) [ pe.reads1 ; pe.reads2 ]
          in
          let tag = id_concat [s.id ; s.species] in
          (*Build biopython index*)
          let index_concat_fasta = Misc_workflows.build_biopythonindex ~tag concat_fasta in
          (*build overlapping read cluster*)
          let cluster_repo = Cdhitoverlap.cdhitoverlap ~tag concat_fasta in
          let rep_cluster_fasta = Cdhitoverlap.cluster_rep cluster_repo in
          let cluster = Cdhitoverlap.cluster cluster_repo in
          (*reformat cluster*)
          let reformated_cluster = reformat_cdhit_cluster ~tag cluster in
          (*build index for cluster*)
          let index_cluster = Misc_workflows.build_biopythonindex ~tag reformated_cluster in
          (*Build blast db of cluster representatives*)
          let parse_seqids = true in
          let hash_index = true in
          let dbtype = "nucl" in
          let cluster_rep_blast_db = BlastPlus.makeblastdb ~hash_index ~parse_seqids ~dbtype  (s.id ^ "_" ^ s.species) rep_cluster_fasta in
          Some {
            Apytram.s ; concat_fasta; index_concat_fasta;
            rep_cluster_fasta; reformated_cluster; index_cluster ;
            cluster_rep_blast_db
          }
        else
          None
      )

  module Seq_dispatcher = struct
    let seq_dispatcher
        ?s2s_tab_by_family
        ~ref_db
        ~query
        ~query_species
        ~query_id
        ~ref_transcriptome
        ~threads
        ~seq2fam : [`seq_dispatcher] directory =
      let open Shell_dsl in
      Workflow.shell ~np:threads ~version:9 ~descr:("SeqDispatcher.py:" ^ query_id ^ "_" ^ query_species) [
        mkdir_p tmp;
        cmd "python" ~img:caars_img [
          file_dump (string Scripts.seq_dispatcher);
          option (flag string "--sp2seq_tab_out_by_family" ) s2s_tab_by_family;
          opt "-d" ident (seq ~sep:"," (List.map ref_db ~f:(fun blast_db -> seq [dep blast_db ; string "/db"]) ));
          opt "-tmp" ident tmp ;
          opt "-log" seq [ dest ; string ("/SeqDispatcher." ^ query_id ^ "." ^ query_species ^ ".log" )] ;
          opt "-q" dep query ;
          opt "-qs" string query_species ;
          opt "-qid" string query_id ;
          opt "-threads" ident np ;
          opt "-t" dep ref_transcriptome ;
          opt "-t2f" dep seq2fam;
          opt "-out" seq [ dest ; string ("/Trinity." ^ query_id ^ "." ^ query_species )] ;
        ]
      ]

    let get_fasta dir (sample : Rna_sample.t) family =
      let file = sprintf "Trinity.%s.%s.%s.fa" sample.id sample.species family in
      Workflow.select dir [file]

  end


  let trinity_annotated_families_of_trinity_assemblies config_dir assemblies ref_blast_dbs threads =
    assoc_map assemblies  ~f:(fun (s : Rna_sample.t) trinity_assembly ->
        let ref_db = List.map s.reference_species ~f:(fun r -> ref_blast_dbs $ r) in
        let query = trinity_assembly in
        let query_species= s.species in
        let query_id = s.id in
        let tag = id_concat s.reference_species in
        let ref_transcriptome =
          List.map s.reference_species ~f:(Configuration_directory.ref_transcriptome config_dir)
          |> Misc_workflows.fasta_concat ~tag:(tag ^ ".ref_transcriptome")
        in
        let seq2fam =
          List.map s.reference_species ~f:(Configuration_directory.ref_seq_fam_links config_dir)
          |> Misc_workflows.fasta_concat ~tag:(tag ^ ".seq2fam")
        in
        Seq_dispatcher.seq_dispatcher
          ~s2s_tab_by_family:true
          ~query
          ~query_species
          ~query_id
          ~ref_transcriptome
          ~seq2fam
          ~ref_db
          ~threads
      )



  let concat_without_error ?(descr="") l : fasta file =
    let open Bistro.Shell_dsl in
    let script =
      let vars = [
        "FILE", seq ~sep:"" l ;
        "DEST", dest ;
      ]
      in
      Commons.bash_script vars {|
        touch tmp
        cat tmp $FILE > tmp1
        mv tmp1 $DEST
        |}
    in
    Workflow.shell ~descr:("concat_without_error" ^ descr) [
      mkdir_p tmp;
      cd tmp;
      cmd "sh" [ file_dump script];
    ]

  let build_target_query dataset ref_species family (trinity_annotated_fams : [`seq_dispatcher] directory Rna_sample.assoc) apytram_group =
    let seq_dispatcher_results_dirs =
      assoc_opt (Dataset.apytram_samples dataset) ~f:(fun s ->
          if String.(s.group_id = apytram_group) && Poly.(s.reference_species = ref_species) && s.run_trinity then
            Some (List.Assoc.find_exn ~equal:Poly.( = ) trinity_annotated_fams s)
          else
            None
        )
    in
    let tag = family ^ ".seqdispatcher" in
    List.map seq_dispatcher_results_dirs ~f:(fun (s, dir) ->
        Seq_dispatcher.get_fasta dir s family
      )
    |> Misc_workflows.fasta_concat ~tag


  let apytram_annotated_ref_fams_by_fam_by_groups (dataset : Dataset.t) configuration_dir trinity_annotated_fams reads_blast_dbs memory_per_sample =
    let apytram_groups = Dataset.apytram_groups dataset in
    let apytram_ref_species = Dataset.apytram_reference_species dataset in
    let apytram_samples = Dataset.apytram_samples dataset in
    List.map dataset.used_families ~f:(fun fam ->
        let fws =
          List.concat_map apytram_groups ~f:(fun apytram_group ->
              let pairs = List.cartesian_product apytram_ref_species [fam] in
              List.concat_map pairs ~f:(fun (ref_species, fam) ->
                  let tag = id_concat [fam.name ; String.concat ~sep:"_" ref_species ; String.strip apytram_group] in
                  let guide_query =
                    List.map ref_species ~f:(fun sp -> Configuration_directory.ref_fams configuration_dir sp fam.name)
                    |> Misc_workflows.fasta_concat ~tag
                  in
                  let target_query = build_target_query dataset ref_species fam.name trinity_annotated_fams apytram_group in
                  let query = Misc_workflows.fasta_concat ~tag:(tag ^ ".+seqdispatcher") [guide_query; target_query] in
                  let compressed_reads_dbs = List.filter_map reads_blast_dbs ~f:(fun ((s : Rna_sample.t), db) ->
                      if (List.equal String.equal s.reference_species ref_species && String.equal s.group_id apytram_group)
                      then Some db else None
                    )
                  in
                  let time_max = 18000 * List.length compressed_reads_dbs in
                  let w =
                    Apytram.apytram_multi_species
                      ~descr:tag ~time_max ~no_best_file:true ~write_even_empty:true
                      ~mal:66 ~i:5 ~evalue:1e-10 ~out_by_species:true
                      ~memory:memory_per_sample ~fam:fam.name ~query compressed_reads_dbs in
                  List.filter_map apytram_samples ~f:(fun s ->
                      if List.equal String.equal s.reference_species ref_species && String.(s.group_id = apytram_group) then
                        Some (s, fam, Apytram.get_fasta w ~family_name:fam.name ~sample_id:s.id )
                      else None
                    )
                )
              )
        in
        (fam, fws)
      )

  let make ?(memory = 4) ?(nthreads = 2) (dataset : Dataset.t) =
    let memory_per_sample, threads_per_sample =
      let nb_samples = List.length dataset.samples in
      Int.(max 1 (memory / (max 1 nb_samples))), Stdlib.(max 1 (nthreads / (max 1 nb_samples)))
    in
    (* let memory_per_thread = Int.(max 1 (config.memory / config.nthreads)) in *)
    let config_dir = Configuration_directory.make ~memory dataset in
    let ref_blast_dbs = ref_blast_dbs dataset config_dir in
    let fasta_reads = fasta_reads dataset in
    let normalized_fasta_reads = normalize_fasta_reads fasta_reads memory_per_sample memory threads_per_sample in
    let trinity_assemblies = trinity_assemblies_of_norm_fasta normalized_fasta_reads ~memory:memory_per_sample ~nthreads:threads_per_sample in
    let trinity_orfs = transdecoder_orfs_of_trinity_assemblies trinity_assemblies ~memory:memory_per_sample ~nthreads:threads_per_sample in
    let trinity_assemblies_stats = assemblies_stats_of_assemblies trinity_assemblies in
    let trinity_orfs_stats = assemblies_stats_of_assemblies trinity_orfs in
    let trinity_annotated_fams = trinity_annotated_families_of_trinity_assemblies config_dir trinity_orfs ref_blast_dbs threads_per_sample in
    let reads_blast_dbs = blast_dbs_of_norm_fasta normalized_fasta_reads in
    let apytram_annotated_ref_fams_by_fam_by_groups =
      apytram_annotated_ref_fams_by_fam_by_groups dataset config_dir trinity_annotated_fams reads_blast_dbs memory_per_sample
    in
    { ref_blast_dbs ; fasta_reads ; normalized_fasta_reads ;
      trinity_assemblies ; trinity_orfs ; trinity_assemblies_stats ;
      trinity_orfs_stats ; trinity_annotated_fams ;
      reads_blast_dbs ; apytram_annotated_ref_fams_by_fam_by_groups }
end
