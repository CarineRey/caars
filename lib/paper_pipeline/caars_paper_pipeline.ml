open Core
open Bistro
(*open Biotope*)
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

let minReadLen_ids = function
  | `Mouse -> 50
  | `Stickleback -> 100

let img_sratoolkit = [ docker_image ~account:"pegi3s" ~name:"sratoolkit" ~tag:"2.10.0" () ]
let img_caars_prep_data = [ docker_image ~account:"carinerey" ~name:"caars_tuto_prep_data" ~tag:"latest" () ]

let sra_of_input = function
    | `id id -> string id
    | `idw w -> string_dep w
    | `file w -> dep w

let fastq_dump_pe_gz ~minReadLen input =
    let sra = sra_of_input input in
    let dir =
      Workflow.shell ~descr:"sratoolkit.fastq_dump" [
        mkdir_p dest ;
        cmd ~img:img_sratoolkit "fastq-dump" [
          opt "-O" Fn.id dest ;
          string "--split-files" ;
          string "--gzip" ;
          opt "--defline-seq" string "'@$ac_$si/$ri'";
          opt "--defline-qual" string "'+'";
          opt "-N" int 10000 ; (*TODO add a debug tag to not download the whole sra*)
          opt "-X" int 200000 ;
          opt "-M" int minReadLen ;
          sra ;
        ] ;
        mv (dest // "*_1.fastq*") (dest // "reads_1.fq.gz") ;
        mv (dest // "*_2.fastq*") (dest // "reads_2.fq.gz") ;
      ]
    in
    Workflow.select dir ["reads_1.fq.gz"],
    Workflow.select dir ["reads_2.fq.gz"]

let fastq_gz srrid =
  fastq_dump_pe_gz ~minReadLen:50 (`id srrid) (*TODO change minReadLen in fonction of the species (maybe useless, if trivial yes else no)*)

let concat_fastq_gz fqs : fastq file =
  let open Bistro.Shell_dsl in
  Workflow.shell ~descr:"concat_fastq_gz" [
    cmd "cat" ~stdout:dest [
      list ~sep:" " Bistro_unix.Cmd.gzdep fqs
    ]
  ]

let pair_map (x, y) ~f = (f x, f y)

let concatenated_fq species =
  srr_ids species
  |> List.map ~f:fastq_gz
  |> List.unzip
  |> pair_map ~f:concat_fastq_gz

let mouse_fq_1, mouse_fq_2 = concatenated_fq `Mouse
let stickleback_fq_1, stickleback_fq_2 = concatenated_fq `Stickleback

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
  
let input_data_repo = Bistro_utils.Repo.[
    item ["sample_sheet.tsv"] sample_sheet ;
    item ["tree.nw"] tree ;
    item ["msa"] (msas_ensemblcompara ()) ;
    item ["rna_seq" ; "Mus_musculus_1.fq"] mouse_fq_1 ;
    item ["rna_seq" ; "Mus_musculus_2.fq"] mouse_fq_2 ;
    item ["rna_seq" ; "Gasterosteus_aculeatus_1.fq"] stickleback_fq_1 ;
    item ["rna_seq" ; "Gasterosteus_aculeatus_2.fq"] stickleback_fq_2 ;
  ]

let main ?(np = 4) ?(memory = 4) ~outdir () =
  let open Bistro_utils in
  let loggers = [
    Console_logger.create () ;
  ]
  in
  Repo.build_main ~loggers ~np ~mem:(`GB memory) ~outdir:(Filename.concat outdir "input_data") input_data_repo

let command =
  let open Command.Let_syntax in
  Command.basic
    ~summary:"caars_paper_pipeline"
    [%map_open
      let outdir = flag "--outdir" (required string) ~doc:"PATH Destination directory."
      and np = flag "--np" (optional int) ~doc:"INT Number of CPUs "
      and memory = flag "--memory" (optional int)    ~doc:"INT Number of GB of system memory to use.(Default:1)"
      in
      main ?np ?memory ~outdir
    ]
