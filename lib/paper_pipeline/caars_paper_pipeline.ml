open Core
open Bistro
open Biotope
open Biotope.Formats
open Bistro.Shell_dsl

(* paired-end data from Mus musculus kidney [(Fushan AA et al., "Gene
   expression defines natural changes in mammalian lifespan.", Aging
   Cell, 2015 Feb
   9;14(3):352-65)](http://www.ncbi.nlm.nih.gov/pubmed/25677554) *)
let mus_musculus_kidney_sample_srr_ids = ["SRR636916";"SRR636917";"SRR636918"]

(* paired-end data from Gasterosteus aculeatus kidney (Stickleback) *)
let stickleback_srr_ids = ["SRR528539" ; "SRR528540"]

let srr_ids = function
  | `Mouse -> mus_musculus_kidney_sample_srr_ids
  | `Stickleback -> stickleback_srr_ids

let minReadLen = function
  | `Mouse -> 50
  | `Stickleback -> 100

let img_caars_prep_data = [ docker_image ~account:"carinerey" ~name:"caars_tuto_prep_data" ~tag:"latest" () ]

let fastq_gz ?(preview = false) species srrid =
  let _N_, _X_ =  if preview then Some 10_000, Some 200_000 else None, None in
  Sra_toolkit.(
    fastq_dump_pe fastq_gz
      ?_N_ ?_X_
      ~defline_qual:"+"
      ~defline_seq:"@$ac_$si/$ri"
      ~minReadLen:(minReadLen species)
      (`id srrid)
  )

let concat_fastq_gz fqs : fastq file =
  let open Bistro.Shell_dsl in
  Workflow.shell ~descr:"concat_fastq_gz" [
    cmd "cat" ~stdout:dest [
      list ~sep:" " Bistro_unix.Cmd.gzdep fqs
    ]
  ]

let pair_map (x, y) ~f = (f x, f y)

let concatenated_fq ?preview species =
  srr_ids species
  |> List.map ~f:(fastq_gz ?preview species)
  |> List.unzip
  |> pair_map ~f:concat_fastq_gz

let get_msas_ensemblcompara_script : text file =
  Bistro_unix.wget "https://github.com/CarineRey/caars/wiki/src/Get_MSAs_EnsemblCompara.pl"

let msas_ensemblcompara ?(preview = false) () : [`msa_ensemblcompara] directory =
  let families = if preview then "ENSGT00550000074800,ENSGT00550000074846,ENSGT00870000136549" else "all" in
  let open Shell_dsl in
  let quoted_string s = quote ~using:'"' (string s) in
  Workflow.shell ~descr:"get_msas_ensemblcompara" [
    cmd "perl" ~img:img_caars_prep_data [
      dep get_msas_ensemblcompara_script ;
      opt "-s" quoted_string "Anolis_carolinensis,Taeniopygia_guttata,Dasypus_novemcinctus,Loxodonta_africana,Sus_scrofa,Ovis_aries,Myotis_lucifugus,Felis_catus,Mustela_putorius_furo,Otolemur_garnettii,Callithrix_jacchus,Homo_sapiens,Oryctolagus_cuniculus,Cavia_porcellus,Ictidomys_tridecemlineatus,Danio_rerio,Lepisosteus_oculatus" ;
      opt "-f" string families ;
      opt "-g" quoted_string "Homo_sapiens" ;
      opt "-r" quoted_string "Mus_musculus,Cavia_porcellus,Gasterosteus_aculeatus,Danio_rerio,Homo_sapiens,Ictidomys_tridecemlineatus,Oryzias_latipes,Oreochromis_niloticus" ;
      opt "-o" Fn.id dest ;
    ]
  ]

let tree : newick file =
  let open Shell_dsl in
  Workflow.shell ~descr:"tree" [
    cmd "cat" ~stdout:dest [
      file_dump (string "(((Danio_rerio,Gasterosteus_aculeatus),Lepisosteus_oculatus),((Anolis_carolinensis,Taeniopygia_guttata),((Loxodonta_africana,Dasypus_novemcinctus),(((((Ictidomys_tridecemlineatus,Mus_musculus),Cavia_porcellus),Oryctolagus_cuniculus),(Otolemur_garnettii,(Callithrix_jacchus,Homo_sapiens))),((Myotis_lucifugus,(Felis_catus,Mustela_putorius_furo)),(Sus_scrofa,Ovis_aries))))));\n") ;
    ]
  ]

let sample_sheet_string =
"id\tspecies\tgroup\tref_species\tpath_fastq_single\tpath_fastq_left\tpath_fastq_right\torientation\trun_trinity\tpath_assembly\trun_apytram\n
CMM\tMus_musculus\tg1\tHomo_sapiens\t-\tinput_data/rna_seq/Mus_musculus_1.fq\tinput_data/rna_seq/Mus_musculus_2.fq\tUP\tyes\t-\tyes\n
CGA\tGasterosteus_aculeatus\tg2\tHomo_sapiens\t-\tinput_data/rna_seq/Gasterosteus_aculeatus_1.fq\tinput_data/rna_seq/Gasterosteus_aculeatus_2.fq\tUP\tyes\t-\tyes\n
"

let sample_sheet : text file =
  let open Shell_dsl in
  Workflow.shell ~descr:"sample_sheeet" [
    cmd "cat" ~stdout:dest [
      file_dump (string sample_sheet_string) ;
    ]
  ]
  
let input_data_repo ?preview () =
  let mouse_fq_1, mouse_fq_2 = concatenated_fq ?preview `Mouse in
  let stickleback_fq_1, stickleback_fq_2 = concatenated_fq ?preview `Stickleback in
  Bistro_utils.Repo.[
    item ["sample_sheet.tsv"] sample_sheet ;
    item ["tree.nw"] tree ;
    item ["msa"] (msas_ensemblcompara ?preview ()) ;
    item ["rna_seq" ; "Mus_musculus_1.fq"] mouse_fq_1 ;
    item ["rna_seq" ; "Mus_musculus_2.fq"] mouse_fq_2 ;
    item ["rna_seq" ; "Gasterosteus_aculeatus_1.fq"] stickleback_fq_1 ;
    item ["rna_seq" ; "Gasterosteus_aculeatus_2.fq"] stickleback_fq_2 ;
  ]

let prepare_dataset ?(np = 4) ?(memory = 4) ?preview ~outdir () =
  let open Bistro_utils in
  let loggers = [
    Console_logger.create () ;
  ]
  in
  let repo = input_data_repo ?preview () in
  Repo.build_main ~loggers ~np ~mem:(`GB memory) ~outdir:(Filename.concat outdir "input_data") repo

let prepare_dataset_command =
  let open Command.Let_syntax in
  Command.basic
    ~summary:"Dataset preparation"
    [%map_open
      let outdir = flag "--outdir" (required string) ~doc:"PATH Destination directory."
      and np = flag "--np" (optional int) ~doc:"INT Number of CPUs "
      and memory = flag "--memory" (optional int)    ~doc:"INT Number of GB of system memory to use.(Default:1)"
      and preview = flag "--preview" no_arg ~doc:" Preview mode"
      in
      prepare_dataset ?np ?memory ~outdir ~preview
    ]

let analysis ?(np = 4) ?(memory = 4) ?family_to_use ?(ali_sister_threshold = 0.0) ~indir ~outdir ~just_parse_input () =
  let indir fn = Filename.concat indir fn in
  let sample_sheet = indir "sample_sheet.tsv" in
  let species_tree_file = indir "tree.nw" in
  let species_tre_filename_abs = if Filename.is_relative species_tree_file then (Filename.concat (Sys.getcwd ()) species_tree_file) else species_tree_file in
  let alignments_dir = indir "msa/MSA/" in
  let seq2sp_dir = indir "msa/Seq2SpTable/" in
  let no_reconcile = false in
  let ali_sister_threshold = ali_sister_threshold in
  let merge_criterion = "merge" in
  let debug = false in
  let get_reads = false in
  let just_parse_input = just_parse_input in
  let html_report = outdir ^ "/report.html" in
  let quiet = false in
  let use_docker = true in
  let family_to_use = family_to_use in
  Caars.App.main
    ~sample_sheet ~species_tree_file:species_tre_filename_abs ~alignments_dir ~ali_sister_threshold ~seq2sp_dir
    ~no_reconcile ~merge_criterion ~debug ~get_reads ~just_parse_input ~html_report
    ~quiet ~use_docker ?family_to_use ~np ~memory ~outdir
    ()

let analysis_command =
  let open Command.Let_syntax in
  Command.basic
    ~summary:"Caars pipeline"
    [%map_open
      let outdir = flag "--outdir" (required string) ~doc:"PATH Destination directory"
      and indir = flag "--indir" (required string) ~doc:"PATH Input directory"
      and np = flag "--np" (optional int) ~doc:"INT Number of CPUs "
      and memory = flag "--memory" (optional int) ~doc:"INT Number of GB of system memory to use.(Default:1)"
      and just_parse_input = flag "--just-parse-input"  no_arg ~doc:" Parse input and exit. Recommended to check all input files. (Default:false)"
      and family_to_use = flag "--family-subset"  (optional Filename.arg_type)    ~doc:"PATH A file containing a subset of families to use.  Default: off"
      and ali_sister_threshold = flag "--mpast"           (optional float)  ~doc:"FLOAT Minimal percentage of alignment of a caars sequence on its (non Caars) closest sequence to be kept in the final output"
      in
      analysis ?np ?memory ?family_to_use ?ali_sister_threshold ~outdir ~indir ~just_parse_input
    ]
