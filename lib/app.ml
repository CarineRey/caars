open Core
open Bistro

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

let main sample_sheet outdir species_tree_file alignments_dir seq2sp_dir np memory no_reconcile _refinetree (*refineali*) ali_sister_threshold merge_criterion debug (get_reads:bool) just_parse_input html_report quiet use_docker family_to_use () =
  let open Defs in
  let open Bistro_utils in
  let loggers quiet html_report = [
    if quiet then Bistro_engine.Logger.null else Console_logger.create () ;
    (match html_report with
     | Some path -> Html_logger.create path
     | None -> Bistro_engine.Logger.null) ;
  ]
  in
  let run_reconciliation = not no_reconcile in
  let ali_sister_threshold = Option.value ~default:0.0 ali_sister_threshold in
  let merge_criterion =
    Option.value_map ~default:Merge merge_criterion ~f:(fun s ->
        match merge_criterion_of_string s with
        | Some mgc -> mgc
        | None -> failwith "Invalid merge criterion argument"
      )
  in
  let nthreads = Option.value ~default:2 np in
  let memory = Option.value ~default:1 memory in
  let dataset = Dataset.make ?family_subset_file:family_to_use ~sample_sheet ~species_tree_file ~alignments_dir ~seq2sp_dir () in
  let pipeline = Pipeline.make ~nthreads ~memory ~merge_criterion dataset ~filter_threshold:ali_sister_threshold ~refine_ali:false ~run_reconciliation in
  let caars_workflow =
    if just_parse_input then
      just_parse_workflow ~outdir pipeline
    else
      full_analysis_workflow ~get_reads ~debug ~outdir pipeline
  in
  Bistro_engine.Scheduler.simple_eval_exn
    ~loggers:(loggers quiet html_report)
    ~np:nthreads
    ~mem:(`GB memory)
    ~collect:true
    ~db_path:"_caars"
    ~allowed_containers:(if use_docker then [`Docker] else [])
    caars_workflow

let command =
  let open Command.Let_syntax in
  Command.basic
    ~summary:"caars"
    [%map_open
      let sample_sheet = flag "--sample-sheet"    (required Filename.arg_type)   ~doc:"PATH sample sheet file."
      and outdir = flag "--outdir"          (required string) ~doc:"PATH Destination directory."
      and species_tree_file = flag "--species-tree"    (required Filename.arg_type)   ~doc:"ABSOLUTE PATH Species tree file in nw format containing all species. Warning absolute path is required."
      and alignments_dir = flag "--alignment-dir"   (required string) ~doc:"PATH Directory containing all gene family alignments (Family_name.fa) in fasta format."
      and seq2sp_dir = flag "--seq2sp-dir"      (required string) ~doc:"PATH Directory containing all linked files (Family_name.tsv). Each line contains a sequence name and its species spaced by a tabulation."
      and np = flag "--np"              (optional int)    ~doc:"INT Number of CPUs (at least 2). (Default:2)"
      and memory = flag "--memory"          (optional int)    ~doc:"INT Number of GB of system memory to use.(Default:1)"
      and no_reconcile = flag "--no-reconcile"    no_arg            ~doc:" Not run final Reconciliation step"
      and refinetree = flag "--refinetree"      no_arg            ~doc:" Refine topology during final Reconciliation step (Default:false)"
      (*  and = flag "--refineali"       no_arg            ~doc:"Refine MSA after the final Reconciliation step (Default:false)"*)
      and ali_sister_threshold = flag "--mpast"           (optional float)  ~doc:"FLOAT Minimal percentage of alignment of a caars sequence on its (non Caars) closest sequence to be kept in the final output"
      and merge_criterion = flag "--merge-criterion" (optional string) ~doc:"STR Merge criterion during redundancy removing. It must be “length“ or “length_complete” or “merge”. “length” means the longest sequence is selected. “length.complete” : means the largest number of complete sites (no gaps). “merge” means that the set of monophyletic sequences is used to build one long “chimera” sequence corresponding to the merging of them."
      and debug = flag "--debug"           no_arg            ~doc:" Get intermediary files (Default:false)"
      and get_reads = flag "--get-reads"       no_arg            ~doc:" Get normalized reads (Default:false)"
      and just_parse_input = flag "--just-parse-input"no_arg            ~doc:" Parse input and exit. Recommended to check all input files. (Default:false)"
      and html_report = flag "--html-report"    (optional string)  ~doc:"PATH Logs build events in an HTML report"
      (* and = flag "--dag-graph"      (optional string)  ~doc:"PATH Write dag graph in a dot file (Can take a lot of time)" *)
      and quiet = flag "--quiet"           no_arg            ~doc:" Do not report progress.  Default: off"
      and use_docker = flag "--use-docker"      no_arg            ~doc:" Use docker in caars.  Default: off"
      and family_to_use = flag "--family-subset"  (optional Filename.arg_type)    ~doc:"PATH A file containing a subset of families to use.  Default: off"
      in
      main sample_sheet outdir species_tree_file alignments_dir seq2sp_dir np memory no_reconcile refinetree (*refineali*) ali_sister_threshold merge_criterion debug get_reads just_parse_input html_report quiet use_docker family_to_use
    ]
