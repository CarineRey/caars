open Core_kernel
open Defs

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
