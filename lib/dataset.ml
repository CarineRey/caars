open Core_kernel
open Defs

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
