type t = {
  id : string ;
  group_id : string ;
  species : string ;
  reference_species : string list ;
  sample_file : Sample_source.t ;
  run_trinity : bool ;
  run_transdecoder : bool ;
  run_apytram : bool ;
  precomputed_assembly : string option ;
}
type 'a assoc = (t * 'a) list
