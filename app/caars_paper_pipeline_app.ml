open Core

let () =
  Command.group ~summary:"CAARS's paper pipeline" [
    "prepare-dataset", Caars_paper_pipeline.prepare_dataset_command ;
    "analysis", Caars_paper_pipeline.analysis_command ;
  ]
  |> Command.run
