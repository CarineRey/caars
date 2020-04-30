open Core_kernel
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
