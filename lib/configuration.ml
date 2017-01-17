open Core.Std
open Commons

type rna_sample = {
  id : string ;
  species : string ;
  ref_species : string ;
  sample_fastq : string sample_fastq ;
  run_trinity : bool ;
  run_transdecoder : bool ;
  path_assembly : string ;
  given_assembly : bool ;
  run_apytram : bool ;
}

type config_rna_seq = rna_sample list

type t = {
  config_rna_seq : config_rna_seq ;
  apytram_samples: rna_sample list ;
  trinity_samples : rna_sample list ;
  all_ref_samples : rna_sample list ;
  all_ref_species : string list ;
  all_apytram_ref_species : string list ;
  families : string list;
  sample_sheet : string ;
  species_tree_file : string ;
  alignments_dir : string ;
  seq2sp_dir : string ;
  outdir : string ;
  threads : int;
  memory : int;
}

let parse_fastq_path = function
  | "-" -> None
  | x -> Some x

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
  | [ id ; species ; ref_species ; path_fastq_single ; path_fastq_left ; path_fastq_right ; orientation ; run_trinity ; path_assembly ; run_apytram] ->
     let run_transdecoder = true in
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

     let sample_fastq = match (parse_fastq_path path_fastq_single,
                             parse_fastq_path path_fastq_left,
                             parse_fastq_path path_fastq_right,
                             parse_orientation id orientation,
                             run_apytram,
                             run_trinity) with
       | ( None    , Some  _, Some _, Some (Right o), _    , _    ) -> Fastq_Paired_end (path_fastq_left, path_fastq_right, o)
       | ( Some  _ , None   , None  , Some (Left o) , _    , _    ) -> Fastq_Single_end (path_fastq_single, o)
       | ( None    , None   , None  , None          , false, true ) -> Fastq_Single_end ("-", US)
       | ( None    , None   , None  , _             , true , _    ) -> failwith ({|Syntax error in the sample sheet file (sample -> |} ^ id ^ {|): You didn't give any RNA-seq data, but you ask to run apytram : it is impossible, apytram needs raw RNA-seq data.|})
       | ( Some  _ , None   , None  , Some (Right o), _    , _    ) -> failwith ({|Syntax error in the sample sheet file (sample -> |} ^ id ^ {|): Incompatible choice. Path for a single-end data but an orientation for a paired-end data.|})
       | ( None    , Some  _, Some _, Some (Left o) , _    , _    ) -> failwith ({|Syntax error in the sample sheet file (sample -> |} ^ id ^ {|): Incompatible choice. Paths for paired-end data but an orientation for a single-end data.|})
       | ( _       , _      , _     , None          , _    , _    ) -> failwith ({|Syntax error in the sample sheet file (sample -> |} ^ id ^ {|): No given orientation.|})
       | _ -> failwith ({|Syntax error in the sample sheet file (sample -> |} ^ id ^ {|): Incompatible choices. path_fastq_single must be "-" if data are "paired-end" and path_fastq_left and path_fastq_right must be "-" if data are "single-end".|})(*(path_fastq_single ^ path_fastq_left ^ path_fastq_right ^ orientation)*)
     in
     { id ;
       species ;
       ref_species ;
       sample_fastq ;
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
    if Filename.check_suffix f ".fa" || Filename.check_suffix f ".fasta" then
      true
    else
      (printf "Warning: %s is not a fasta file (extention must be .fa or .fasta)\n" f ; false)
    )
  |> Array.map ~f:(fun f -> fst (String.lsplit2_exn f ~on:'.')) (* Il y a obligatoirement un point dans le nom du fichier fasta *)
  |> Array.to_list


let load ~sample_sheet ~species_tree_file ~alignments_dir ~seq2sp_dir ~np ~memory ~outdir =
  let threads = match np with
    | x when x > 1 -> np
    | _ -> failwith "The number of CPUs must be at least 2"
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
        Some s.ref_species
      else
        None
    )
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

  if List.contains_dup id_list then
    failwith {|There are duplicate id in the first colum of the config file.|}
  else if Filename.is_relative species_tree_file then
    failwith {|amalgam needs the absolute path of the species tree.|}
  else
    {
      config_rna_seq;
      apytram_samples = filter_samples (fun s -> s.run_apytram);
      trinity_samples = filter_samples (fun s -> s.run_trinity);
      all_ref_samples = filter_samples (fun s -> s.run_apytram || s.run_trinity);
      all_ref_species = filter_ref_species;
      all_apytram_ref_species = filter_apytram_ref_species;
      families = families_of_alignments_dir alignments_dir;
      sample_sheet ;
      species_tree_file ;
      alignments_dir ;
      seq2sp_dir ;
      threads;
      memory ;
      outdir;
    }
