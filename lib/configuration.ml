open Core
open Commons

type t = {
  config_rna_seq : config_rna_seq ;
  apytram_samples: rna_sample list ;
  trinity_samples : rna_sample list ;
  all_ref_samples : rna_sample list ;
  all_ref_species : string list ;
  all_apytram_ref_species : string list list;
  families : string list;
  sample_sheet : string ;
  species_tree_file : string ;
  alignments_dir : string ;
  seq2sp_dir : string ;
  outdir : string ;
  threads : int;
  memory : int;
  run_reconciliation : bool;
  refinetree : bool;
  refineali : bool;
  debug : bool;
  just_parse_input : bool;
  ali_sister_threshold : float;
}

let parse_fastq_path = function
  | "-" -> None
  | x -> Some x

let parse_fastX_path f = match ( f, Filename.split_extension f) with
  | ("-", _ ) -> None
  | (x, ( _, Some "fa"  ))-> Some (Fasta x)
  | (x, ( _, Some "fasta"  )) -> Some (Fasta x)
  | (x, ( _, Some "fq"  )) -> Some (Fastq x)
  | (x, ( _, Some "fastq"  )) -> Some (Fastq x)
  | (x, ( _, Some y  )) -> failwith ({|Syntax error: sample file extension must be in  ["fa", ".fasta","fq","fastq"] (detected extension: |} ^ y  ^ ",  " ^ x ^" )." )
  |  _  -> failwith ({|Syntax error: sample file extension must be in  ["fa", "fasta","fq","fastq"] |})

let parse_orientation id = function
  | "F"  -> Some (Left F)
  | "R"  -> Some (Left R)
  | "US" -> Some (Left US)
  | "RF" -> Some (Right RF)
  | "FR" -> Some (Right FR)
  | "UP" -> Some (Right UP)
  | "-"  -> None
  | _    -> failwith ({|Syntax error in the sample sheet file (sample -> |} ^ id ^ {|: orientation must be in ["F","R","RF","FR","US","UP"] |})

let parse_line_fields_of_rna_conf_file = function
  | [ id ; species ; ref_species ; path_fastx_single ; path_fastx_left ; path_fastx_right ; orientation ; run_trinity ; path_assembly ; run_apytram] ->
     let run_transdecoder = true in

     let ref_species = List.sort compare (String.split ~on:',' ref_species) in
     let run_trinity = match run_trinity with
       | "yes" | "Yes" | "y" | "Y" -> true
       | "no" | "No" | "n" | "N" -> false
       | _ -> failwith ({| Syntax error inthe sample sheet file (sample -> |} ^ id ^ {|): run_trinity must be "yes" or "no" |})
     in
     let (path_assembly,given_assembly) = match (path_assembly,run_trinity) with
       | ("-" ,_) -> ("-",false)
       | (path,true) -> (path,true)
       | (path,false) -> failwith ({|Syntax error in the sample sheet file (sample -> |} ^ id ^ {|: you gave a trinity assembly path but run_trinity is false. It is incompatible.|})
     in
     let run_apytram = match run_apytram with
       | "yes" | "Yes" | "y" | "Y" -> true
       | "no" | "No" | "n" | "N" -> false
       | _ -> failwith ({|Syntax error in the sample sheet file (sample -> |} ^ id ^ {|: run_apytram must be "yes" or "no" |})
     in
     let sample_file= match (parse_fastX_path path_fastx_single,
                             parse_fastX_path path_fastx_left,
                             parse_fastX_path path_fastx_right,
                             parse_orientation id orientation,
                             run_apytram,
                             run_trinity) with
       | ( None    , Some  (Fastq _ ), Some (Fastq _ ), Some (Right o), _    , _   ) -> Sample_fastq (Fastq_Paired_end (path_fastx_left, path_fastx_right, o))
       | ( None    , Some  (Fasta _ ), Some (Fasta _ ), Some (Right o), _    , _   ) -> Sample_fasta (Fasta_Paired_end (path_fastx_left, path_fastx_right, o))
       | ( Some  (Fastq _) , None   , None  , Some (Left o) , _    , _    ) -> Sample_fastq(Fastq_Single_end (path_fastx_single, o))
       | ( Some  (Fasta _) , None   , None  , Some (Left o) , _    , _    ) -> Sample_fasta(Fasta_Single_end (path_fastx_single, o))
       | ( None    , None   , None  , None  , false, true ) -> Sample_fastq (Fastq_Single_end ("-", US))
       | ( None    , Some  (Fasta _ ), Some (Fastq _ ), _ , _  , _   ) -> failwith ({|Syntax error in the sample sheet file (sample -> |} ^ id ^ {|): You give  a fasta and a fastq file |})
       | ( None    , Some  (Fastq _ ), Some (Fasta _ ), _ , _  , _   ) -> failwith ({|Syntax error in the sample sheet file (sample -> |} ^ id ^ {|): You give  a fasta and a fastq file |})
       | ( None    , None   , None  , _             , true , _    ) -> failwith ({|Syntax error in the sample sheet file (sample -> |} ^ id ^ {|): You didn't give any RNA-seq data, but you ask to run apytram : it is impossible, apytram needs raw RNA-seq data.|})
       | ( Some  _ , None   , None  , Some (Right o), _    , _    ) -> failwith ({|Syntax error in the sample sheet file (sample -> |} ^ id ^ {|): Incompatible choice. Path for a single-end data but an orientation for a paired-end data.|})
       | ( None    , Some  _, Some _, Some (Left o) , _    , _    ) -> failwith ({|Syntax error in the sample sheet file (sample -> |} ^ id ^ {|): Incompatible choice. Paths for paired-end data but an orientation for a single-end data.|})
       | ( _       , _      , _     , None          , _    , _    ) -> failwith ({|Syntax error in the sample sheet file (sample -> |} ^ id ^ {|): No given orientation.|})
       | _ -> failwith ({|Syntax error in the sample sheet file (sample -> |} ^ id ^ {|): Incompatible choices. path_fastq_single must be "-" if data are "paired-end" and path_fastq_left and path_fastq_right must be "-" if data are "single-end".|})(*(path_fastq_single ^ path_fastq_left ^ path_fastq_right ^ orientation)*)
     in
     { id ;
       species ;
       ref_species ;
       sample_file ;
       run_trinity ;
       run_transdecoder ;
       path_assembly ;
       given_assembly ;
       run_apytram
     }
  | _ -> failwith "Syntax error in the sample sheet file. There aren't 10 tab delimited columns."

let parse_rna_conf_file path =
  In_channel.read_lines path
  |> List.tl_exn
  |> List.map ~f:(String.split ~on:'\t')
  |> List.map ~f:parse_line_fields_of_rna_conf_file


let families_of_alignments_dir alignments_dir =
  Sys.readdir alignments_dir
  |> Array.filter ~f:(fun f ->
    if Filename.check_suffix f ".fa" then
      true
    else
      (printf "Warning: %s is not a fasta file (extention must be .fa)\n" f ; false)
    )
  |> Array.map ~f:(fun f -> fst (String.lsplit2_exn f ~on:'.')) (* Il y a obligatoirement un point dans le nom du fichier fasta *)
  |> Array.to_list


let load ~sample_sheet ~species_tree_file ~alignments_dir ~seq2sp_dir ~np ~memory ~run_reconciliation ~refinetree ~refineali ~ali_sister_threshold ~debug ~just_parse_input ~outdir =
  let threads = match (np, run_reconciliation) with
    | (x, true) when x > 1 -> np
    | (x, false) when x > 0 -> np
    | (x, true) -> failwith "The number of CPUs must be at least 2 if you want to use reconciliation"
    | (x, false) -> failwith "The number of CPUs must be at least 1"
  in
  let memory = match np with
    | x when x > 0 -> memory
    | _ -> failwith "The memory must be at least 1"
  in
  let config_rna_seq = parse_rna_conf_file sample_sheet in
  let filter_samples f = List.filter config_rna_seq ~f in
  let id_list = List.map config_rna_seq ~f:(fun s -> s.id) in
  let filter_ref_species = List.filter_map config_rna_seq ~f:(fun s ->
      if (s.run_apytram || s.run_trinity) then
        Some (s.ref_species:string list)
      else
        None
    )
    |> List.concat
  in
  let filter_apytram_ref_species =
    let all_sorted = List.sort compare (List.filter_map config_rna_seq ~f:(fun s ->
        if s.run_apytram then
          Some s.ref_species
        else
          None
      )) in
    let rec uniq l = match l with
      | x :: y :: z when x = y -> uniq (x :: z)
      | x :: y :: z -> x :: uniq (y ::z)
      | [x] -> [x]
      | [] -> []
    in
    uniq all_sorted
  in
  let families = families_of_alignments_dir alignments_dir in
  let _ = (printf "%i families.\n" (List.length families); ())  in

  if List.contains_dup id_list then
    failwith {|There are duplicate id in the first colum of the config file.|}
  else if Filename.is_relative species_tree_file then
    failwith {|caars needs the absolute path of the species tree.|}
  else if families = [] then
    failwith ({|No files with .fa extention in |} ^ alignments_dir)
  else
    {
      config_rna_seq;
      apytram_samples = filter_samples (fun s -> s.run_apytram);
      trinity_samples = filter_samples (fun s -> s.run_trinity);
      all_ref_samples = filter_samples (fun s -> s.run_apytram || s.run_trinity);
      all_ref_species = filter_ref_species;
      all_apytram_ref_species = filter_apytram_ref_species;
      families;
      sample_sheet ;
      species_tree_file ;
      alignments_dir ;
      seq2sp_dir ;
      threads;
      memory ;
      run_reconciliation ;
      refinetree ;
      refineali ;
      debug ;
      just_parse_input ;
      outdir ;
      ali_sister_threshold ;
    }
