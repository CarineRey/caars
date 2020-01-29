open Core
open Commons

type t = {
  config_rna_seq : config_rna_seq ;
  apytram_samples: rna_sample list ;
  trinity_samples : rna_sample list ;
  all_ref_samples : rna_sample list ;
  all_ref_species : string list ;
  all_apytram_ref_species : string list list;
  apytram_group_list : string list ;
  all_families : family list;
  used_families : family list;
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
  get_reads : bool;
  just_parse_input : bool;
  ali_sister_threshold : float;
  merge_criterion : merge_criterion;
}


(* let parse_fastq_path = function
 *   | "-" -> None
 *   | x -> Some x *)

let parse_fastX_path f = match ( f, Filename.split_extension f) with
  | ("-", _ ) -> None
  | (x, ( _, Some "fa"  ))-> Some (Fasta x)
  | (x, ( _, Some "fasta"  )) -> Some (Fasta x)
  | (x, ( _, Some "fq"  )) -> Some (Fastq x)
  | (x, ( _, Some "fastq"  )) -> Some (Fastq x)
  | (x, ( _, Some y  )) -> failwith ({|Syntax error: sample file extension must be in  ["fa", ".fasta","fq","fastq"] (detected extension: |} ^ y  ^ ",  " ^ x ^" )." )
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
    | h::t -> "col #" ^ (string_of_int  i) ^ ": " ^ h ^ ";" ^ (str_elements (i+1) t)
  in
  "[" ^ (str_elements 1 l) ^ "]"

let parse_line_fields_of_rna_conf_file = function
  | [ id ; species ; apytram_group ; ref_species ; path_fastx_single ; path_fastx_left ; path_fastx_right ; orientation ; run_trinity ; path_assembly ; run_apytram] ->
     let run_transdecoder = true in
     let ref_species = List.sort ~compare:String.compare (String.split ~on:',' ref_species) in
     let run_trinity = match run_trinity with
       | "yes" | "Yes" | "y" | "Y" -> true
       | "no" | "No" | "n" | "N" -> false
       | _ -> failwith ({| Syntax error inthe sample sheet file (sample -> |} ^ id ^ {|): run_trinity must be "yes" or "no" |})
     in
     let (path_assembly,given_assembly) = match (path_assembly,run_trinity) with
       | ("-", _) -> ("-",false)
       | (path, true) -> (path,true)
       | (_, false) -> failwith ({|Syntax error in the sample sheet file (sample -> |} ^ id ^ {|: you gave a trinity assembly path but run_trinity is false. It is incompatible.|})
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
       | ( Some  _ , None   , None  , Some (Right _), _    , _    ) -> failwith ({|Syntax error in the sample sheet file (sample -> |} ^ id ^ {|): Incompatible choice. Path for a single-end data but an orientation for a paired-end data.|})
       | ( None    , Some  _, Some _, Some (Left _) , _    , _    ) -> failwith ({|Syntax error in the sample sheet file (sample -> |} ^ id ^ {|): Incompatible choice. Paths for paired-end data but an orientation for a single-end data.|})
       | ( _       , _      , _     , None          , _    , _    ) -> failwith ({|Syntax error in the sample sheet file (sample -> |} ^ id ^ {|): No given orientation.|})
       | _ -> failwith ({|Syntax error in the sample sheet file (sample -> |} ^ id ^ {|): Incompatible choices. path_fastq_single must be "-" if data are "paired-end" and path_fastq_left and path_fastq_right must be "-" if data are "single-end".|})(*(path_fastq_single ^ path_fastq_left ^ path_fastq_right ^ orientation)*)
     in
     Some { id ;
       species ;
       apytram_group ;
       ref_species ;
       sample_file ;
       run_trinity ;
       run_transdecoder ;
       path_assembly ;
       given_assembly ;
       run_apytram
     }
  | [""] -> (printf "Warning: empty line in the sample sheet file\n"; None)
  | l -> failwith ("Syntax error in the sample sheet file. There aren't 11 tab delimited columns: " ^ (str_list_sample_line l))

module SS = Set.Make(String)

let not_include l1 l2 =
    let rec include_rec l1 l2 att it = match (l1, l2, att) with
        | h1 :: t1, h2 :: t2, att when Poly.(h1 = h2) -> include_rec t1 t2 att it
        | [],  _ , _-> it
        | h1 :: t1, h2 :: t2, att -> include_rec (h1 :: t1)  t2 (h2 :: att) it
        | h1 :: t1, [] , att -> include_rec t1  att [] (h1 :: it)
    in
    include_rec l1 l2 [] []

let print_list l =
  let rec str_elements = function
    | [] -> ""
    | h::t -> h ^ "; " ^ (str_elements t)
  in
  (str_elements l)

let parse_family_to_use_file all_families path =
   let fs = In_channel.read_lines path
    |> List.map ~f:(String.strip)
    |> List.filter_map ~f:(function
        | "" -> None
        | x -> Some x)
    in
    let sorted_fs = List.sort ~compare:String.compare fs in
    let sorted_all_families = List.sort ~compare:String.compare all_families in
    let diff = not_include sorted_fs sorted_all_families in
    match diff with
        | [] -> sorted_fs
        | l-> failwith ("Families [" ^ (print_list l) ^ "] don't exist (or not unique). See in " ^ path)


let parse_rna_conf_file path =
  In_channel.read_lines path
  |> List.tl_exn (* remove the first line*)
  |> List.map ~f:(String.split ~on:'\t')
  |> List.filter_map ~f:parse_line_fields_of_rna_conf_file

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


let load ~sample_sheet ~species_tree_file ~alignments_dir ~seq2sp_dir ~np ~memory ~run_reconciliation ~refinetree ~refineali ~ali_sister_threshold ~merge_criterion ~debug ~get_reads ~just_parse_input ~outdir ~family_to_use =
    let threads = match (np, run_reconciliation) with
    | (x, true) when x > 1 -> np
    | (x, false) when x > 0 -> np
    | (_, true) -> failwith "The number of CPUs must be at least 2 if you want to use reconciliation"
    | (_, false) -> failwith "The number of CPUs must be at least 1"
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
  let filter_apytram_ref_species, apytram_group_list =
    let rec uniq l = match l with
      | x :: y :: z when Poly.(x = y) -> uniq (x :: z)
      | x :: y :: z -> x :: uniq (y ::z)
      | [x] -> [x]
      | [] -> []
    in
    let all_sorted =
      List.filter_map config_rna_seq ~f:(fun s ->
          if s.run_apytram then
            Some s.ref_species
          else
            None
        )
      |> List.sort ~compare:Poly.compare in
    let all_group =
      List.filter_map config_rna_seq ~f:(fun s ->
          if s.run_apytram then
            Some s.apytram_group
          else
            None
        )
      |> List.sort ~compare:Poly.compare in

    uniq all_sorted, uniq all_group
  in
  let all_families_noid = families_of_alignments_dir alignments_dir in

  let used_families_noid = match family_to_use with
    | None -> all_families_noid
    | Some path -> parse_family_to_use_file all_families_noid path

  in

  let atribute_id fam_l =
    let rec att_id fam_l f_id res = match fam_l with
      | name :: t -> att_id t (f_id+1) ({name; f_id} :: res)
      | [] -> res
    in
    att_id fam_l 1 []
  in
  let all_families = atribute_id all_families_noid in

  let used_families = List.map used_families_noid ~f:(fun u_f ->
      let id = List.filter_map  all_families ~f:(fun fam ->
        if String.(fam.name = u_f) then
          Some fam.f_id
        else
          None)
      |> List.hd
      in
      {name = u_f; f_id = (match id with
                          | Some x -> x
                          | _ ->  0 );}
  )
  in
  let _ = (printf "%i families in %s.\n" (List.length all_families) alignments_dir; ())  in
  let _ = (printf "%i families will be used.\n" (List.length used_families); ())  in
  let merge_criterion = parse_merge_criterion merge_criterion in

  let _ = if debug then (printf "debug: %b.\n" (debug))  else ()  in
  let _ = if get_reads then (printf "get_reads: %b.\n" (get_reads))  else ()  in

  if List.contains_dup id_list ~compare:String.compare then
    failwith {|There are duplicate id in the first colum of the config file.|}
  else if Filename.is_relative species_tree_file then
    failwith {|caars needs the absolute path of the species tree.|}
  else if Poly.(all_families = []) then
    failwith ({|No files with .fa extention in |} ^ alignments_dir)
  else
    {
      config_rna_seq;
      apytram_samples = filter_samples (fun s -> s.run_apytram);
      trinity_samples = filter_samples (fun s -> s.run_trinity);
      all_ref_samples = filter_samples (fun s -> s.run_apytram || s.run_trinity);
      all_ref_species = filter_ref_species;
      all_apytram_ref_species = filter_apytram_ref_species;
      apytram_group_list ;
      all_families;
      used_families;
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
      get_reads ;
      just_parse_input ;
      outdir ;
      ali_sister_threshold ;
      merge_criterion ;
    }
