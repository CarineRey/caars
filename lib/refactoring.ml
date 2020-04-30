open Core_kernel
open Bistro
open Defs
open Wutils

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
      |> List.filter ~f:String.(( <> ) "")
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

  let trinity_samples dataset =
    List.filter dataset.samples ~f:(fun s -> s.run_trinity)

  let reference_samples dataset =
    List.filter dataset.samples ~f:(fun s -> s.run_apytram || s.run_trinity)

  let has_at_least_one_sample_with_reference dataset =
    not (List.is_empty (reference_samples dataset))
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

  let family_metadata dir = Workflow.select dir ["FamilyMetadata.txt"]
  let species_metadata dir = Workflow.select dir ["SpeciesMetadata.txt"]
  let usable_families dir = Workflow.select dir ["UsableFamilies.txt"]
  let detected_families dir = Workflow.select dir ["DetectedFamilies.txt"]
end

module Pipeline = struct
  type t = {
    run_reconciliation : bool ;
    configuration_directory : Configuration_directory.t ;
    checked_used_families_all_together : text file ;
    fasta_reads : fasta file OSE_or_PE.t Rna_sample.assoc ;
    normalized_fasta_reads : fasta file OSE_or_PE.t Rna_sample.assoc ;
    trinity_assemblies : fasta file Rna_sample.assoc ;
    trinity_orfs : fasta file Rna_sample.assoc ;
    trinity_assemblies_stats : text file Rna_sample.assoc ;
    trinity_orfs_stats : text file Rna_sample.assoc ;
    trinity_annotated_fams : [`seq_dispatcher] directory Rna_sample.assoc ;
    ref_blast_dbs : blast_db file assoc ;
    reads_blast_dbs : Apytram.compressed_read_db Rna_sample.assoc ;
    apytram_orfs_ref_fams : (Family.t * (Rna_sample.t * Family.t * fasta file) list) list ;
    apytram_checked_families : (Family.t * (Rna_sample.t * Family.t * fasta file) list) list ;
    apytram_annotated_families : (Family.t * fasta file) list ;
    merged_families : (Family.t * [ `seq_integrator ] directory * [ `seq_integrator ] directory option) list ;
    merged_and_reconciled_families : (Family.t * Generax.phylotree directory * [ `seq_integrator ] directory) list ;
    merged_reconciled_and_realigned_families_dirs : [`merged_families_distributor] directory ;
    reconstructed_sequences : [`reconstructed_sequences] directory option ;
    orthologs_per_seq : [`extract_orthologs] directory ;
    final_plots : [`final_plots] directory ;
  }

  let rna_sample_needs_rna (s : Rna_sample.t) =
    match s.run_apytram, s.run_trinity, s.precomputed_assembly with
    | true, _, _          -> true
    | false, true, Some _ -> false
    | false, true, None   -> true
    | false, false, _     -> false

  let check_used_families ~used_fam_list ~usable_fam_file =
    let open Bistro.Shell_dsl in
    let sorted_usable_fam_file = tmp // "usablefam.sorted.txt" in
    let sorted_all_used_fam_file = dest in
    let common_fam_file = tmp // "common_fam.txt" in
    let fam_subset_not_ok = tmp // "fam_subset_not_ok.txt" in
    let all_used_fam = Bistro.Template_dsl.(
        List.map used_fam_list ~f:(fun (fam : Family.t) -> seq [string fam.name])
        |> seq ~sep:"\n"
      )
    in

    let script_post ~fam_subset_not_ok =
      let args = [
        "FILE_EMPTY", fam_subset_not_ok ;
      ]
      in
      Commons.bash_script args {|
    if [ -s $FILE_EMPTY ]
    then
      echo "These families are not in the \"Usable\" families:"
      cat $FILE_EMPTY
      echo "Use the option --just-parse-input and --family-subset with an empty file to get the file UsableFamilies.txt"
      exit 3
    else
      exit 0
    fi
    |}
    in
    Workflow.shell ~descr:("check_used_families") [
      mkdir_p tmp;
      cmd "sort" ~stdout:sorted_usable_fam_file [ dep usable_fam_file;];
      cmd "sort" ~stdout:sorted_all_used_fam_file [ file_dump all_used_fam; ];
      cmd "join" ~stdout: common_fam_file [ string "-1 1"; sorted_all_used_fam_file; sorted_usable_fam_file];
      cmd "comm" ~stdout: fam_subset_not_ok [string "-3"; common_fam_file;sorted_all_used_fam_file];
      cmd "bash" [ file_dump (script_post ~fam_subset_not_ok)]
    ]

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

    let fasta_file_name (sample : Rna_sample.t) family =
      sprintf "Trinity.%s.%s.%s.fa" sample.id sample.species family

    let get_fasta dir (sample : Rna_sample.t) family =
      Workflow.select dir [fasta_file_name sample family]
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

  (* This is needed by [build_target_query] to concat a list of fasta
     in a directory without error, even if some requested files are
     not present. It seems that seq_dispatcher doesn't produce files
     for all families, hence the need to do this. *)
  let concat_without_error ?tag l : fasta file =
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
    Workflow.shell ~descr:(descr ?tag "concat_without_error") [
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
    concat_without_error ~tag (
      List.map seq_dispatcher_results_dirs ~f:Bistro.Shell_dsl.(fun (s, dir) ->
          dep dir // Seq_dispatcher.fasta_file_name s family
        )
    )

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

  let checkfamily
      ?(descr="")
      ~ref_db
      ~(input:fasta file)
      ~family
      ~ref_transcriptome
      ~seq2fam
      ~evalue
    : fasta file =
    let open Bistro.Shell_dsl in
    let tmp_checkfamily = tmp // "tmp" in
    let dest_checkfamily = dest // "sequences.fa" in

    Workflow.shell ~version:8 ~descr:("CheckFamily.py" ^ descr) [
      mkdir_p tmp_checkfamily;
      cd tmp_checkfamily;
      cmd "python" ~img:caars_img [
        file_dump (string Scripts.check_family);
        opt "-tmp" ident tmp_checkfamily ;
        opt "-i" dep input;
        opt "-t" dep ref_transcriptome ;
        opt "-f" string family;
        opt "-t2f" dep seq2fam;
        opt "-o" ident dest_checkfamily;
        opt "-d" ident (seq ~sep:"," (List.map ref_db ~f:(fun blast_db -> seq [dep blast_db ; string "/db"]) ));
        opt "-e" float evalue;
      ]
    ]
    |> Fn.flip Workflow.select [ "sequences.fa" ]

  let apytram_checked_families_of_orfs_ref_fams apytram_orfs_ref_fams configuration_dir ref_blast_dbs =
    List.map apytram_orfs_ref_fams ~f:(fun (fam, fws) ->
        let checked_fws = List.map fws ~f:(fun ((s : Rna_sample.t), (f : Family.t), apytram_orfs_fasta) ->
            let input = apytram_orfs_fasta in
            let tag = id_concat s.reference_species in
            let ref_transcriptome =
              List.map s.reference_species ~f:(Configuration_directory.ref_transcriptome configuration_dir)
              |> Misc_workflows.fasta_concat ~tag:(tag ^ ".ref_transcriptome")
            in
            let seq2fam =
              List.map s.reference_species ~f:(Configuration_directory.ref_seq_fam_links configuration_dir)
              |> Misc_workflows.fasta_concat ~tag:(tag ^ ".seq2fam") in
            let ref_db =
              List.map s.reference_species ~f:(( $ ) ref_blast_dbs) in
            let checked_families_fasta =
              checkfamily ~descr:(":"^s.id^"."^f.name) ~input ~family:f.name ~ref_transcriptome ~seq2fam ~ref_db ~evalue:1e-6
            in
            (s, f, checked_families_fasta)
          ) in
        (fam, checked_fws)
      )

  let parse_apytram_results apytram_annotated_ref_fams =
    let open Bistro.Shell_dsl in
    List.map apytram_annotated_ref_fams ~f:(fun (fam, fws) ->
        let config = Bistro.Template_dsl.(
            List.map fws ~f:(fun ((s : Rna_sample.t), (f : Family.t), w) ->
                seq ~sep:"\t" [ string s.species ; string s.id ; string f.name ; int f.id ; dep w ]
              )
            |> seq ~sep:"\n"
          )
        in
        let fw =
          Workflow.shell ~version:4 ~descr:("parse_apytram_results.py." ^ fam.Family.name) ~np:1  [
            cmd "python" ~img:caars_img [
              file_dump (string Scripts.parse_apytram_results) ;
              file_dump config ;
              dest ]
          ]
        in
        (fam, fw)
      )

  module Seq_integrator = struct
    open Bistro.Shell_dsl

    let transform_species_list l =
      list ~sep:"," string l

    let seq_integrator
        ?realign_ali
        ?resolve_polytomy
        ?species_to_refine_list
        ?no_merge
        ?merge_criterion
        ~family
        ~trinity_fam_results_dirs
        ~apytram_results_dir
        ~alignment_sp2seq
        alignment
      : [`seq_integrator] directory =

      let open Bistro.Shell_dsl in
      let merge_criterion_string =
        Option.bind merge_criterion ~f:(function
            | Merge -> None
            | Length -> Some "length"
            | Length_complete -> Some "length.complete"
          )
      in
      let get_trinity_file_list extension dirs =
        List.concat_map dirs ~f:(fun ((s : Rna_sample.t), dir) ->
            [ dep dir ; string ("/Trinity." ^ s.id ^ "." ^ s.species ^ "." ^ family ^ "." ^ extension) ; string ","]
          )
      in
      let get_apytram_file_list extension dir =
        [ dep dir ; string ("/apytram." ^ family ^ "." ^ extension) ; string ","]
      in
      let trinity_fasta_list = get_trinity_file_list "fa" trinity_fam_results_dirs in
      let trinity_sp2seq_list = get_trinity_file_list "sp2seq.txt" trinity_fam_results_dirs in

      let apytram_fasta = get_apytram_file_list "fa" apytram_results_dir in
      let apytram_sp2seq = get_apytram_file_list "sp2seq.txt" apytram_results_dir in

      let sp2seq = List.concat [[dep alignment_sp2seq ; string "," ] ; trinity_sp2seq_list ; apytram_sp2seq ]  in
      let fasta = List.concat [trinity_fasta_list; apytram_fasta]  in

      let tmp_merge = tmp // "tmp" in

      Workflow.shell ~version:12 ~descr:("SeqIntegrator.py:" ^ family) [
        mkdir_p tmp_merge ;
        cmd "python" ~img:caars_img [
          file_dump (string Scripts.seq_integrator);
          opt "-tmp" ident tmp_merge;
          opt "-log" seq [ tmp_merge ; string ("/SeqIntegrator." ^ family ^ ".log" )] ;
          opt "-ali" dep alignment ;
          opt "-fa" (seq ~sep:"") fasta;
          option (flag string "--realign_ali") realign_ali;
          option (opt "--merge_criterion" string) merge_criterion_string;
          option (flag string "--no_merge") no_merge;
          option (flag string "--resolve_polytomy") resolve_polytomy;
          opt "-sp2seq" (seq ~sep:"") sp2seq  ; (* list de sp2seq delimited by comas *)
          opt "-out" seq [ dest ; string "/" ; string family] ;
          option (opt "-sptorefine" transform_species_list) species_to_refine_list;
        ]
      ]

    let seq_filter
        ?realign_ali
        ?resolve_polytomy
        ?species_to_refine_list
        ~filter_threshold
        ~family
        ~alignment
        ~tree
        ~sp2seq
      : [`seq_integrator] directory  =

      let open Bistro.Shell_dsl in
      let tmp_merge = tmp // "tmp" in

      Workflow.shell ~version:8 ~descr:("SeqFilter.py:" ^ family) [
        mkdir_p tmp_merge ;
        cmd "python" ~img:caars_img [
          file_dump (string Scripts.seq_filter);
          opt "-tmp" ident tmp_merge;
          opt "-log" seq [ tmp_merge ; string ("/SeqFilter." ^ family ^ ".log" )] ;
          opt "-ali" dep alignment ;
          opt "-t" dep tree;
          opt "--filter_threshold" float filter_threshold;
          option (flag string "--realign_ali") realign_ali;
          option (flag string "--resolve_polytomy") resolve_polytomy;
          opt "-sp2seq" dep sp2seq  ;
          opt "-out" seq [ dest ; string "/" ; string family] ;
          option (opt "-sptorefine" transform_species_list) species_to_refine_list;
        ]
      ]

    let alignment dir (fam : Family.t) = Workflow.select dir [ fam.name ^ ".fa" ]
    let tree dir (fam : Family.t) = Workflow.select dir [ fam.name ^ ".tree" ]
    let sp2seq dir (fam : Family.t) = Workflow.select dir [ fam.name ^ ".sp2seq.txt" ]

  end


  let merged_families_of_families (dataset : Dataset.t) configuration_dir trinity_annotated_fams apytram_annotated_fams merge_criterion filter_threshold =
    List.map dataset.used_families ~f:(fun family ->
        let trinity_fam_results_dirs=
          List.map (Dataset.trinity_samples dataset) ~f:(fun s ->
              (s , List.Assoc.find_exn ~equal:Poly.( = ) trinity_annotated_fams s)
            )
        in
        let apytram_results_dir =  List.Assoc.find_exn ~equal:Poly.( = ) apytram_annotated_fams family in
        let alignment = Workflow.input (dataset.alignments_dir ^ "/" ^ family.name ^ ".fa")  in
        let alignment_sp2seq = Configuration_directory.ali_species2seq_links configuration_dir family.name  in
        let species_to_refine_list = List.map (Dataset.reference_samples dataset) ~f:(fun s -> s.species) in
        let w = if (List.length species_to_refine_list) = 0 then
            Seq_integrator.seq_integrator ~realign_ali:false ~resolve_polytomy:true ~no_merge:true ~family:family.name ~trinity_fam_results_dirs ~apytram_results_dir ~alignment_sp2seq ~merge_criterion alignment
          else
            Seq_integrator.seq_integrator ~realign_ali:false ~resolve_polytomy:true ~species_to_refine_list ~family:family.name ~trinity_fam_results_dirs ~apytram_results_dir ~alignment_sp2seq ~merge_criterion alignment
        in
        let tree = Seq_integrator.tree w family in
        let alignment = Seq_integrator.alignment w family in
        let sp2seq = Seq_integrator.sp2seq w family in

        let wf =
          if List.length species_to_refine_list > 0 then
            Some (Seq_integrator.seq_filter ~realign_ali:true ~resolve_polytomy:true ~filter_threshold ~species_to_refine_list ~family:family.name ~tree ~alignment ~sp2seq)
          else None
        in
        (family, w, wf )
      )

  let generax_by_fam_of_merged_families (dataset : Dataset.t) merged_families memory threads =
    List.map  merged_families ~f:(fun ((fam : Family.t), merged_without_filter_family, merged_and_filtered_family) ->
        let merged_family = match merged_and_filtered_family with
          | Some w -> w
          | None -> merged_without_filter_family
        in

        let tree = Seq_integrator.tree merged_family fam in
        let alignment = Seq_integrator.alignment merged_family fam in
        let sp2seq = Seq_integrator.sp2seq merged_family fam in
        let sptreefile = Workflow.input dataset.species_tree_file in
        (fam, Generax.generax ~family:fam.name ~descr:(":" ^ fam.name) ~threads ~memory ~sptreefile ~link:sp2seq ~tree alignment, merged_family)
      )

  let merged_families_distributor dataset merged_reconciled_and_realigned_families ~run_reconciliation ~refine_ali : [`merged_families_distributor] directory =
    let open Bistro.Shell_dsl in
    let more_than_one_sample_with_reference = Dataset.has_at_least_one_sample_with_reference dataset in
    let extension_list_merged = [(".fa","out/MSA_out");(".tree","out/GeneTree_out");(".sp2seq.txt","no_out/Sp2Seq_link")] in
    let extension_list_filtered = [(".discarded.fa","out/FilterSummary_out");(".filter_summary.txt","out/FilterSummary_out")] in

    let extension_list_reconciled = [("_ReconciledTree.nw","","out/GeneTreeReconciled_nw");
                                     ("_ReconciledTree.nhx", "", "out/GeneTreeReconciled_out");
                                     (".events.txt", "", "out/DL_out");
                                     (".orthologs.txt", "", "out/Orthologs_out")] in
    let dest_dir_preparation_commands = List.concat [
        [
          mkdir_p tmp;
          mkdir_p (dest // "out" // "MSA_out");
          mkdir_p (dest // "out" // "GeneTree_out");
          mkdir_p (dest // "no_out" // "Sp2Seq_link");
        ] ;

        if more_than_one_sample_with_reference
        then [ mkdir_p (dest // "out" // "FilterSummary_out") ]
        else [] ;

        if run_reconciliation then
          [
            mkdir_p (dest // "out" // "GeneTreeReconciled_out");
            mkdir_p (dest // "out" // "DL_out");
            mkdir_p (dest // "out" // "Orthologs_out");
          ]
        else [] ;

        if refine_ali && run_reconciliation then
          [mkdir_p (dest // "Realigned_fasta")]
        else []
      ]
    in
    let commands_for_one_family ((f : Family.t), reconciled_w, merged_w) =
      let open Bistro.Template_dsl in
      List.concat [
        List.map extension_list_merged ~f:(fun (ext,dir) ->
            let input = Workflow.select merged_w [ f.name ^ ext ] in
            let output = dest // dir // (f.name ^ ext)  in
            seq ~sep:" " [ string "cp"; dep input ; ident output ]
          ) ;
        if more_than_one_sample_with_reference then
          List.map extension_list_filtered ~f:(fun (ext,dir) ->
              let input = Workflow.select merged_w [ f.name  ^ ext ] in
              let output = dest // dir // (f.name  ^ ext)  in
              seq ~sep:" " [ string "cp"; dep input ; ident output ]
            )
        else [] ;
        if run_reconciliation then
          List.concat [
            List.map extension_list_reconciled ~f:(fun (ext,dirin,dirout) ->
                let input = Workflow.select reconciled_w [ dirin ^ f.name  ^ ext ] in
                let output = dest // dirout // (f.name  ^ ext)  in
                seq ~sep:" " [ string "cp"; dep input ; ident output ]
              )
            ;
          ]
        else [] ;
      ]
    in
    let script =
      List.concat_map merged_reconciled_and_realigned_families ~f:commands_for_one_family
      |> seq ~sep:"\n"
    in
    let commands = dest_dir_preparation_commands @ [ cmd "bash" [ file_dump script ] ] in
    Workflow.shell ~descr:"build_output_directory" ~version:1 commands

  let get_reconstructed_sequences dataset merged_and_reconciled_families_dirs =
    let open Bistro.Shell_dsl in
    if Dataset.has_at_least_one_sample_with_reference dataset then
      let species_to_refine_list = List.map (Dataset.reference_samples dataset) ~f:(fun s -> s.species) in
      Some (Workflow.shell ~descr:"GetReconstructedSequences.py" ~version:6 [
          mkdir_p dest;
          cmd "python" ~img:caars_img [
            file_dump (string Scripts.get_reconstructed_sequences);
            dep merged_and_reconciled_families_dirs // "out/MSA_out";
            dep merged_and_reconciled_families_dirs // "no_out/Sp2Seq_link";
            seq ~sep:"," (List.map species_to_refine_list ~f:(fun sp -> string sp));
            ident dest
          ]
        ])
    else
      None

  let write_orthologs_relationships dataset merged_and_reconciled_families_dirs ~run_reconciliation =
    let ortho_dir,species_to_refine_list =
      if run_reconciliation then
      Some (Workflow.select merged_and_reconciled_families_dirs ["out/Orthologs_out"]),
      Some (List.map (Dataset.reference_samples dataset) ~f:(fun s -> s.species))
      else (None, None)
    in
    let open Bistro.Shell_dsl in
    Workflow.shell ~descr:"ExtractOrthologs.py" ~version:7 [
      mkdir_p dest;
      cmd "python" ~img:caars_img [
        file_dump (string Scripts.extract_orthologs);
        ident dest;
        dep merged_and_reconciled_families_dirs // "no_out/Sp2Seq_link";
        option (opt "" dep) ortho_dir ;
        option (opt "" Seq_integrator.transform_species_list) species_to_refine_list ;
      ]
    ]

  let build_final_plots dataset orthologs_per_seq merged_reconciled_and_realigned_families_dirs ~run_reconciliation =
    let open Bistro.Shell_dsl in
    let formated_target_species =
      match Dataset.reference_samples dataset with
      | [] -> None
      | samples -> Some (
          List.map samples ~f:(fun s ->
              seq ~sep:":" [string s.species ; string s.id]
            )
        )
    in
    let dloutprefix = dest // "D_count" in
    Workflow.shell ~descr:"final_plots.py" ~version:19 (List.concat [
        [mkdir_p dest;
         cmd "python" ~img:caars_img [
           file_dump (string Scripts.final_plots);
           opt "-i_ortho" dep orthologs_per_seq;
           opt "-i_filter" dep (Workflow.select merged_reconciled_and_realigned_families_dirs ["out/"]);
           opt "-o" ident dest;
           option (opt "-t_sp" (seq ~sep:",")) formated_target_species;
         ];
        ];
        if run_reconciliation then
          [cmd "python" ~img:caars_img [
              file_dump (string Scripts.count_dl);
              opt "-o" ident dloutprefix;
              opt "-sp_tree" dep (Workflow.input (dataset.species_tree_file));
              opt "-rec_trees_dir" dep (Workflow.select merged_reconciled_and_realigned_families_dirs ["out/GeneTreeReconciled_out"])
            ];
          ]
        else
          []

      ])

  let make
      ?(memory = 4) ?(nthreads = 2)
      ~merge_criterion ~filter_threshold
      ~refine_ali ~run_reconciliation
      (dataset : Dataset.t) =
    let memory_per_sample, threads_per_sample =
      let nb_samples = List.length dataset.samples in
      Int.(max 1 (memory / (max 1 nb_samples))), Stdlib.(max 1 (nthreads / (max 1 nb_samples)))
    in
    (* let memory_per_thread = Int.(max 1 (config.memory / config.nthreads)) in *)
    let config_dir = Configuration_directory.make ~memory dataset in
    let checked_used_families_all_together =
      check_used_families ~used_fam_list:dataset.used_families ~usable_fam_file:(Configuration_directory.usable_families config_dir)
    in
    let ref_blast_dbs = ref_blast_dbs dataset config_dir in
    let fasta_reads = fasta_reads dataset in
    let normalized_fasta_reads = normalize_fasta_reads fasta_reads memory_per_sample memory threads_per_sample in
    let trinity_assemblies = trinity_assemblies_of_norm_fasta normalized_fasta_reads ~memory:memory_per_sample ~nthreads:threads_per_sample in
    let trinity_orfs = transdecoder_orfs_of_trinity_assemblies trinity_assemblies ~memory:memory_per_sample ~nthreads:threads_per_sample in
    let trinity_assemblies_stats = assemblies_stats_of_assemblies trinity_assemblies in
    let trinity_orfs_stats = assemblies_stats_of_assemblies trinity_orfs in
    let trinity_annotated_fams = trinity_annotated_families_of_trinity_assemblies config_dir trinity_orfs ref_blast_dbs threads_per_sample in
    let reads_blast_dbs = blast_dbs_of_norm_fasta normalized_fasta_reads in
    let apytram_orfs_ref_fams =
      apytram_annotated_ref_fams_by_fam_by_groups dataset config_dir trinity_annotated_fams reads_blast_dbs memory_per_sample
    in
    let apytram_checked_families = apytram_checked_families_of_orfs_ref_fams apytram_orfs_ref_fams config_dir ref_blast_dbs in
    let apytram_annotated_families = parse_apytram_results apytram_checked_families in
    let merged_families = merged_families_of_families dataset config_dir trinity_annotated_fams apytram_annotated_families merge_criterion filter_threshold in
    let merged_and_reconciled_families = generax_by_fam_of_merged_families dataset merged_families memory nthreads in
    let merged_reconciled_and_realigned_families_dirs =
      merged_families_distributor dataset merged_and_reconciled_families ~refine_ali ~run_reconciliation
    in
    let reconstructed_sequences = get_reconstructed_sequences dataset merged_reconciled_and_realigned_families_dirs in
    let orthologs_per_seq = write_orthologs_relationships dataset merged_reconciled_and_realigned_families_dirs ~run_reconciliation in
    let final_plots = build_final_plots dataset orthologs_per_seq merged_reconciled_and_realigned_families_dirs ~run_reconciliation in
    { run_reconciliation ;
      configuration_directory = config_dir ; checked_used_families_all_together ;
      ref_blast_dbs ; fasta_reads ; normalized_fasta_reads ;
      trinity_assemblies ; trinity_orfs ; trinity_assemblies_stats ;
      trinity_orfs_stats ; trinity_annotated_fams ;
      reads_blast_dbs ; apytram_orfs_ref_fams ; apytram_checked_families ;
      apytram_annotated_families ; merged_families ;
      merged_and_reconciled_families ; merged_reconciled_and_realigned_families_dirs ;
      reconstructed_sequences ; orthologs_per_seq ; final_plots }
end

module Repo = struct
  open Bistro_utils.Repo

  let ( ++ ) t h = h :: t
  let ( +? ) t maybe_h =
    match maybe_h with
    | None -> t
    | Some h -> h :: t

  let just_parse_repo (p : Pipeline.t) =
    [
      item ["FamilyMetadata.txt"] (Configuration_directory.family_metadata p.configuration_directory) ;
      item ["SpeciesMetadata.txt"] (Configuration_directory.species_metadata p.configuration_directory) ;
      item ["UsableFamilies.txt"] (Configuration_directory.usable_families p.configuration_directory) ;
      item ["DetectedFamilies.txt"] (Configuration_directory.detected_families p.configuration_directory) ;
      item ["UsedFamilies.txt"] p.checked_used_families_all_together
    ]

  let trinity_assemblies (p : Pipeline.t) =
    List.filter_map p.trinity_assemblies ~f:(fun (s, trinity_assembly) ->
        match s.precomputed_assembly with
        | Some _ -> None
        | None ->
          item ["draft_assemblies" ; "raw_assemblies" ; "Draft_assemblies." ^ s.id ^ "_" ^ s.species ^ ".fa"] trinity_assembly
          |> Option.some
      )

  let trinity_orfs (p : Pipeline.t) =
    List.filter_map p.trinity_orfs ~f:(fun (s, trinity_orf) ->
        match s.precomputed_assembly with
        | Some _ -> None
        | None ->
          item ["draft_assemblies" ; "cds" ; "Draft_assemblies.cds." ^ s.id ^ "_" ^ s.species ^ ".fa"] trinity_orf
          |> Option.some
      )

   let target_to_sample_fasta (s : Rna_sample.t) d = function
    | OSE_or_PE.Single_end se -> [ item [ d ; s.id ^ "_" ^ s.species ^ ".fa" ] se.reads ]
    | Paired_end pe -> [
        item [ d ; s.id ^ "_" ^ s.species ^ ".left.fa" ] pe.reads1 ;
        item [ d ; s.id ^ "_" ^ s.species ^ ".right.fa" ] pe.reads2 ;
      ]

  let sample_fasta_repo (p : Pipeline.t) =
    List.concat [
      List.concat_map p.fasta_reads ~f:(fun (s, sample_fasta) -> target_to_sample_fasta s "rna_seq/raw_fasta" sample_fasta) ;
      List.concat_map p.normalized_fasta_reads ~f:(fun (s,norm_fasta) -> target_to_sample_fasta s "rna_seq/norm_fasta" norm_fasta) ;
    ]

  let debug_repo (p : Pipeline.t) =
    List.concat [
      List.map p.trinity_assemblies_stats ~f:(fun (s, trinity_assembly_stats) ->
          item ["debug" ; "trinity_assembly" ; "trinity_assemblies_stats" ; "Trinity_assemblies." ^ s.id ^ "_" ^ s.species ^ ".stats"] trinity_assembly_stats
        ) ;
      List.map p.trinity_orfs_stats ~f:(fun (s, trinity_orfs_stats) ->
          item ["debug" ; "trinity_assembly" ; "trinity_assemblies_stats" ; "Transdecoder_cds." ^ s.id ^ "_" ^ s.species ^ ".stats"] trinity_orfs_stats
        ) ;
      List.map p.trinity_annotated_fams ~f:(fun (s, trinity_annotated_fams) ->
          item ["debug" ; "trinity_blast_annotation" ; "trinity_annotated_fams" ; s.id ^ "_" ^ s.species ^ ".vs." ^ String.concat ~sep:"_" s.reference_species] trinity_annotated_fams
        ) ;
      List.map p.ref_blast_dbs ~f:(fun (ref_species, blast_db) ->
          [ "debug" ; "trinity_blast_annotation" ; "ref_blast_db" ; ref_species ] %> blast_db
        ) ;
      List.map p.reads_blast_dbs ~f:(fun (s,blast_db) ->
          [ "debug" ; "rna_seq" ;"rep_cluster_blast_db" ; s.id ^ "_" ^ s.species ] %> blast_db.cluster_rep_blast_db
        )
      ;
      List.map p.apytram_orfs_ref_fams ~f:(fun (_, fws) ->
          List.map fws ~f:(fun (s, fam, apytram_result) ->
              [ "debug" ; "apytram_assembly" ; "apytram_results_by_ref_by_group_by_fam" ; fam.name ; s.id ^ "_" ^ s.species ^ ".fa" ] %> apytram_result
            )
        ) |> List.concat
      ;
      List.map p.apytram_checked_families ~f:(fun (_, fws) ->
          List.map fws ~f:(fun (s, fam, apytram_result) ->
              [ "debug" ; "apytram_assembly" ; "apytram_checked_families" ; fam.name ; s.id ^ "_" ^ s.species ^ ".fa"] %> apytram_result
            )
        ) |> List.concat
      ;
      List.map p.apytram_annotated_families ~f:(fun (fam, fw) ->
          [["debug" ; "apytram_assembly" ;"apytram_annotated_sequences"; fam.name ] %> fw]
        ) |> List.concat
      ;
      List.concat_map p.merged_families ~f:(fun (fam, merged_family, merged_and_filtered_family) ->
          match (merged_family, merged_and_filtered_family) with
          | (w1, Some w2) ->  [ [ "debug" ; "merged_families" ; fam.name  ] %> w1; [ "debug" ; "merged_filtered_families" ; fam.name  ] %> w2 ]
          | (w1, None) -> [[ "debug" ; "merged_families" ; fam.name  ] %> w1]
        )
    ]

  let full_repo ?(get_reads = false) ?(debug = false) (p : Pipeline.t) =
    List.concat [
      just_parse_repo p ;
      trinity_assemblies p ;
      trinity_orfs p ;
      ([]
       ++ item ["assembly_results_by_fam" ] (Workflow.select p.merged_reconciled_and_realigned_families_dirs ["out/"])
       ++ item ["all_fam.seq2sp.tsv"]       (Workflow.select p.orthologs_per_seq ["all_fam.seq2sp.tsv"])
       ++ item ["plots"] p.final_plots
       +? Option.map p.reconstructed_sequences ~f:(fun w -> item ["assembly_results_only_seq"] (Workflow.select w ["assemblies/"]))
       +? (
         if p.run_reconciliation then
           Some (item ["all_fam.orthologs.tsv"] (Workflow.select p.orthologs_per_seq ["all_fam.orthologs.tsv"]))
         else None
       )
      ) ;
      if get_reads then sample_fasta_repo p else [] ;
      if debug then debug_repo p else [] ;
    ]

  let precious_repo (p : Pipeline.t) =
    let pi = precious_item in
    let normalized_fasta_samples =
      List.concat_map p.normalized_fasta_reads ~f:(function
          | (_, Single_end se) -> [ pi se.reads ]
          | (_, Paired_end pe) -> [ pi pe.reads1 ; pi pe.reads2 ]
        )
    in
    let reads_blast_dbs =
      List.concat_map p.reads_blast_dbs ~f:(fun (_, x) ->
          [pi x.concat_fasta; pi x.index_concat_fasta; pi x.rep_cluster_fasta; pi x.reformated_cluster; pi x.index_cluster ; pi x.cluster_rep_blast_db ]
        )
    in
    let get_merged_families = function
      | (_, w1, Some w2) -> [pi w1; pi w2]
      | (_, w1, None) -> [pi w1]
    in
    List.concat [
      [pi p.configuration_directory];
      normalized_fasta_samples ;
      List.map p.trinity_assemblies ~f:(fun (_, x) -> pi x) ;
      List.map p.trinity_orfs ~f:(fun (_, x) -> pi x) ;
      reads_blast_dbs ;
      List.map p.trinity_annotated_fams ~f:(fun (_, x) -> pi x) ;
      List.concat_map p.merged_families ~f:get_merged_families;
      List.map p.merged_and_reconciled_families ~f:(fun (_, x, _) -> pi x) ;
      List.concat_map p.merged_and_reconciled_families ~f:(fun (_ , w1, w2) -> [pi w1; pi w2;]) ;
      List.concat_map p.apytram_checked_families ~f:(fun (_, fws) -> List.map fws ~f:(fun (_, _, x) -> pi x)) ;
      List.map p.apytram_annotated_families ~f:(fun (_, x) -> pi x) ;
    ]
end


module Report : sig
val generate :
  trinity_assemblies_stats:(Rna_sample.t * string) list ->
  string ->
  unit
end
=
struct
  open Core
  open Tyxml_html

  let k = txt

  let optint = function
    | None -> k"NA"
    | Some i -> k (Int.to_string i)

  let optfloat = function
    | None -> k"NA"
    | Some i -> k (Float.to_string i)

  (* let svg_from_file fn =
   *   let contents = In_channel.read_all fn in
   *   let src = sprintf "data:image/%s;base64,%s" "svg+xml" contents in
   *   Tyxml_html.img ~src ~alt:"" () *)

  let head t =
    head (title (txt t)) [
      link ~rel:[`Stylesheet] ~href:"http://netdna.bootstrapcdn.com/bootstrap/3.0.2/css/bootstrap.min.css" () ;
      link ~rel:[`Stylesheet] ~href:"http://netdna.bootstrapcdn.com/bootstrap/3.0.2/css/bootstrap-theme.min.css" () ;
      script ~a:[a_src "https://code.jquery.com/jquery.js"] (txt "") ;
      script ~a:[a_src "http://netdna.bootstrapcdn.com/bootstrap/3.0.2/js/bootstrap.min.js"] (txt "") ;
    ]

  let trinity_section trinity_assemblies_stats =
    let table_headers = [
      tr [
        th [ k "Species" ] ;
        th [ k "nb_genes" ] ;
        th [ k "nb_transcripts" ] ;
        th [ k "gc"] ;
        th [ k "n50" ] ;
      ]
    ]
    in
    let foreach_sample (sample, assembly_stats) =
      let { Trinity_stats.n50 ; nb_genes; gc; nb_transcripts } = Trinity_stats.parse assembly_stats in
      tr [
        td [ k sample.Rna_sample.species ] ;
        td [ optint nb_genes ] ;
        td [ optint nb_transcripts ] ;
        td [ optfloat gc] ;
        td [ optint n50 ] ;
      ]
    in
    [
      h2 [ k"Trinity assemblies stats" ] ;
      table ~a:[a_class ["table" ; "table-condensed"]] (List.concat [ table_headers; (List.map trinity_assemblies_stats ~f:foreach_sample)])
    ]

  (* http://ocsigen.org/tyxml/4.0.1/manual/intro*)
  let render ~trinity_assemblies_stats =
    let mytitle = "Caars report" in
    let contents = List.concat [
        [
          h1 [txt "A fabulous title"] ;
          txt "This is a fabulous content." ;
        ] ;
        trinity_section trinity_assemblies_stats ;

      ]

    in
    html
      (head mytitle)
      (body [ div ~a:[a_class ["container"]] contents ])



  let save path doc =
    let buf = Buffer.create 253 in
    let formatter = Format.formatter_of_buffer buf in
    Tyxml_html.pp () formatter doc ;
    Out_channel.with_file path ~f:(fun oc ->
        let contents = Buffer.contents buf in
        Out_channel.output_string oc contents
      )

  let generate ~trinity_assemblies_stats dest =
    let doc = render ~trinity_assemblies_stats in
    save dest doc
end

let just_parse_workflow ~outdir p =
  let module R = Bistro_utils.Repo in
  let repo = R.precious_item p.Pipeline.configuration_directory :: Repo.just_parse_repo p in
  R.to_workflow repo ~outdir

let full_analysis_workflow ?get_reads ?debug ~outdir p =
  let module R = Bistro_utils.Repo in
  let repo =
    List.concat [
      Repo.precious_repo p ;
      Repo.full_repo ?get_reads ?debug p ;
    ]
  in
  let repo_workflow = R.to_workflow repo ~outdir in
  let trinity_assemblies_stats =
    List.map p.trinity_assemblies_stats ~f:(fun (s, w) -> Workflow.(both (data s) (path w)))
    |> Workflow.list
  in
  [%workflow
      let () = [%eval repo_workflow] in
      let trinity_assemblies_stats = [%eval trinity_assemblies_stats] in
      Report.generate
        ~trinity_assemblies_stats
        (Filename.concat outdir "report_end.html")
    ]
