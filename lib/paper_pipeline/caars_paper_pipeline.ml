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

let msas_ensemblcompara () : [`msa_ensemblcompara] directory =
  let open Shell_dsl in
  let quoted_string s = quote ~using:'"' (string s) in
  Workflow.shell ~descr:"get_msas_ensemblcompara" [
    cmd "perl" ~img:img_caars_prep_data [
      dep get_msas_ensemblcompara_script ;
      opt "-s" quoted_string "Anolis_carolinensis,Taeniopygia_guttata,Dasypus_novemcinctus,Loxodonta_africana,Sus_scrofa,Ovis_aries,Myotis_lucifugus,Felis_catus,Mustela_putorius_furo,Otolemur_garnettii,Callithrix_jacchus,Homo_sapiens,Oryctolagus_cuniculus,Cavia_porcellus,Ictidomys_tridecemlineatus,Danio_rerio,Lepisosteus_oculatus" ;
      opt "-f" string "ENSGT00550000074800,ENSGT00550000074846,ENSGT00870000136549" ; (*TODO if debug: "ENSGT00550000074800,ENSGT00550000074846,ENSGT00870000136549" else: "all"*)
      opt "-g" quoted_string "Homo_sapiens" ;
      opt "-r" quoted_string "Mus_musculus,Cavia_porcellus,Gasterosteus_aculeatus,Danio_rerio,Homo_sapiens,Ictidomys_tridecemlineatus,Oryzias_latipes,Oreochromis_niloticus" ;
      opt "-o" Fn.id dest ;
    ]
  ]

let tree : newick file =
  let open Shell_dsl in
  Workflow.shell ~descr:"tree" [
    cmd "cat" ~stdout:dest [
      file_dump (string "(((Danio_rerio,Gasterosteus_aculeatus),Lepisosteus_oculatus),((Anolis_carolinensis,Taeniopygia_guttata),((Loxodonta_africana,Dasypus_novemcinctus),(((((Ictidomys_tridecemlineatus,Mus_musculus),Cavia_porcellus),Oryctolagus_cuniculus),(Otolemur_garnettii,(Callithrix_jacchus,Homo_sapiens))),((Myotis_lucifugus,(Felis_catus,Mustela_putorius_furo)),(Sus_scrofa,Ovis_aries))))))") ;
    ]
  ]

let sample_sheet_string =
{|id      species  group  ref_species     path_fastq_single       path_fastq_left path_fastq_right        orientation     run_trinity    path_assembly   run_apytram
CMM     Mus_musculus  g1   Homo_sapiens  -       rna_seq/Mus_musculus_1.fq        rna_seq/Mus_musculus_2.fq        UP    yes      -       yes
CGA     Gasterosteus_aculeatus  g2  Homo_sapiens  -       rna_seq/Gasterosteus_aculeatus_1.fq      rna_seq/Gasterosteus_aculeatus_2.fq      UP    yes      -       yes
|}

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
    item ["msa"] (msas_ensemblcompara ()) ;
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

let analysis ?(np = 4) ?(memory = 4) ~indir ~outdir () =
  let indir fn = Filename.concat indir fn in
  let sample_sheet = indir "sample_sheet.tsv" in
  let species_tree_file = indir "tree.nw" in
  let alignments_dir = indir "msa" in
  let seq2sp_dir = indir (assert false) in
  let no_reconcile = assert false in
  let ali_sister_threshold = assert false in
  let merge_criterion = assert false in
  let debug = assert false in
  let get_reads = false in
  let just_parse_input = false in
  let html_report = assert false in
  let quiet = true in
  let use_docker = true in
  let family_to_use = assert false in
  Caars.App.main
    ~sample_sheet ~species_tree_file ~alignments_dir ~ali_sister_threshold ~seq2sp_dir
    ~no_reconcile ~merge_criterion ~debug ~get_reads ~just_parse_input ~html_report
    ~quiet ~use_docker ~family_to_use ~np ~memory ~outdir
    ()

let analysis_command =
  let open Command.Let_syntax in
  Command.basic
    ~summary:"Caars pipeline"
    [%map_open
      let outdir = flag "--outdir" (required string) ~doc:"PATH Destination directory"
      and indir = flag "--indir" (required string) ~doc:"PATH Input directory"
      and np = flag "--np" (optional int) ~doc:"INT Number of CPUs "
      and memory = flag "--memory" (optional int)    ~doc:"INT Number of GB of system memory to use.(Default:1)"
      in
      analysis ?np ?memory ~outdir ~indir
    ]
