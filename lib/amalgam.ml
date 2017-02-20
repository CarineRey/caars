open Core.Std
open Bistro.Std
open Bistro.EDSL
open Bistro_bioinfo.Std
open Commons
open Configuration

type output = [ `amalgam_output ]

type tabular

let alignement_fasta fam : (output, fasta) selector =
  selector [ "Alignements" ; fam ^ ".fa" ]

let gene_tree fam : (output, [`newick]) selector =
  selector [ "Gene_trees" ; fam ^ ".tree" ]

type sp2seq_link

let sp2seq_link fam : (output, sp2seq_link) selector =
  selector [ "Sp2Seq_link" ; fam ^ ".sp2seq.txt" ]


type configuration_dir = [ `configuration ]

let parse_input { sample_sheet ; species_tree_file ; alignments_dir ; seq2sp_dir; memory } : configuration_dir directory workflow=
  workflow ~np:1 ~descr:"Parse input" ~version:9 ~mem:(memory * 1024) [
    cmd "ParseInput.py"  [ string sample_sheet ;
                           string species_tree_file;
                           string alignments_dir;
                           string seq2sp_dir;
                           ident dest ;
                         ]
  ]

let ref_transcriptomes species : (configuration_dir, fasta) selector =
  selector ["R_Sp_transcriptomes" ;  species ^ "_transcriptome.fa" ]

let ref_seq_fam_links species : (configuration_dir, tabular) selector =
  selector ["R_Sp_Seq_Fam_links";  species ^ "_Fam_Seq.tsv"  ]

let ref_fams species family =
  selector ["R_Sp_Gene_Families"; species ^ "." ^ family ^ ".fa"]

let ali_species2seq_links family =
  selector ["Alignments_Species2Sequences" ; "alignments." ^  family ^ ".sp2seq.txt" ]

let ref_blast_dbs_of_configuration_dir {all_ref_species} configuration_dir =
  List.map all_ref_species ~f:(fun ref_species ->
    let fasta = configuration_dir / ref_transcriptomes ref_species in
    let parse_seqids = true in
    let dbtype = "nucl" in
    (ref_species, BlastPlus.makeblastdb ~parse_seqids ~dbtype  ("DB_" ^ ref_species) fasta)
    )


let fastq_to_fasta_conversion {all_ref_samples} dep_input =
  List.filter_map all_ref_samples ~f:(fun s ->
      let run_conversion = match (s.run_apytram,s.run_trinity, s.given_assembly) with
        |(true,_,_)         -> true
        |(false,true,true)  -> false
        |(false,true,false) -> true
        |(false,false,_)    -> false
      in
      if run_conversion then
        let sample_fastq = sample_fastq_map input s.sample_fastq in
        let sample_fastq_to_sample_fasta = function
          | Fastq_Single_end (w, o ) -> Fasta_Single_end ( Trinity.fastool  ~dep_input w , o )
          | Fastq_Paired_end (lw, rw , o) -> Fasta_Paired_end ( Trinity.fastool ~dep_input lw , Trinity.fastool ~dep_input rw , o)
        in
        let sample_fasta = sample_fastq_to_sample_fasta sample_fastq in
        Some (s,sample_fasta)
      else
        None
    )

let normalize_fasta fasta_reads {threads ; memory } =
  List.map fasta_reads ~f:(fun (s,fasta_sample) ->
      let max_cov = 50 in
      let normalization_dir = Trinity.fasta_read_normalization max_cov ~threads ~memory fasta_sample in

      let norm_fasta_sample_to_normalization_dir normalization_dir = function
        | Fasta_Single_end (w, o ) -> Fasta_Single_end ( normalization_dir / selector ["single.norm.fa"] , o )
        | Fasta_Paired_end (lw, rw , o) -> Fasta_Paired_end ( normalization_dir / selector ["left.norm.fa"] , normalization_dir / selector ["right.norm.fa"], o )
      in
      (s, norm_fasta_sample_to_normalization_dir normalization_dir fasta_sample )
    )


let trinity_assemblies_of_norm_fasta norm_fasta { memory ; threads ; trinity_samples; apytram_samples} =
  List.concat [
    List.filter_map norm_fasta ~f:(fun (s, norm_fasta) ->
        let threads = if (List.length apytram_samples > 0) && (threads > 1) then
                        threads - 1
                      else
                        threads
                      in
        let memory = if (List.length apytram_samples > 0) && (memory > 4) then
                        Pervasives.(memory * 80 / 100)
                      else
                        memory
                      in
        match (s.run_trinity, s.given_assembly) with
        | (true,false) -> Some (s, Trinity.trinity_fasta ~descr:s.species ~no_normalization:true ~full_cleanup:true ~memory ~threads norm_fasta)
        | (_, _)   -> None
      );
    List.filter_map trinity_samples ~f:(fun s ->
        if s.given_assembly then
          Some (s, input s.path_assembly)
        else
          None
      )
  ]

let transdecoder_orfs_of_trinity_assemblies trinity_assemblies { memory ; threads } =
  List.map trinity_assemblies ~f:(fun (s,trinity_assembly) ->
      match (s.run_transdecoder,s.given_assembly) with
      | (true,false) -> let pep_min_length = 50 in
        let retain_long_orfs = 150 in
        (s, Transdecoder.transdecoder ~descr:("Assembly." ^ s.species) ~retain_long_orfs ~pep_min_length ~only_best_orf:false ~memory ~threads trinity_assembly)
      | (false, _ ) -> (s, trinity_assembly)
      | (true, true) -> (s, trinity_assembly)
    )


let assemblies_stats_of_fasta =
  List.map  ~f:(fun (s,assembly) ->
  (s, Trinity.assembly_stats assembly)
  )


let concat = function
  | [] -> raise (Invalid_argument "fastX concat: empty list")
  | x :: [] -> x
  | fXs ->
    workflow ~descr:"fastX.concat" [
      cmd "cat" ~stdout:dest [ list dep ~sep:" " fXs ]
    ]


let blast_dbs_of_norm_fasta norm_fasta =
  List.filter_map norm_fasta ~f:(fun (s, norm_fasta) ->
      if s.run_apytram then
        let parse_seqids = true in
        let dbtype = "nucl" in
        let fasta_to_norm_fasta_sample = function
          | Fasta_Single_end (w, _ ) -> w
          | Fasta_Paired_end (lw, rw , _) -> concat [ lw ; rw ]
        in
        let fasta = fasta_to_norm_fasta_sample norm_fasta in
        Some (s, BlastPlus.makeblastdb ~parse_seqids ~dbtype  (s.id ^ "_" ^ s.species) fasta)
      else
        None
    )


let seq_dispatcher
    ?s2s_tab_by_family
    ?ref_db
    ~query
    ~query_species
    ~query_id
    ~ref_transcriptome
    ~threads
    ~seq2fam : fasta workflow =
  workflow ~np:threads ~version:9 ~descr:("SeqDispatcher.py:" ^ query_species ^ " ") [
    mkdir_p tmp;
    cmd "SeqDispatcher.py"  [
      option (flag string "--sp2seq_tab_out_by_family" ) s2s_tab_by_family;
      option (opt "-d" (fun blast_db -> seq [dep blast_db ; string "/db"])) ref_db;
      opt "-tmp" ident tmp ;
      opt "-log" seq [ dest ; string ("/SeqDispatcher." ^ query_id ^ "." ^ query_species ^ ".log" )] ;
      opt "-q" dep query ;
      opt "-qs" string query_species ;
      opt "-qid" string query_id ;
      opt "-threads" ident np ;
      opt "-t" dep ref_transcriptome ;
      opt "-t2f" dep seq2fam;
      opt "-out" seq [ dest ; string ("/Trinity." ^ query_id ^ "." ^ query_species )] ;
    ]
  ]

let trinity_annotated_fams_of_trinity_assemblies configuration_dir ref_blast_dbs configuration =
  let threads = configuration.threads in
  List.map ~f:(fun (s,trinity_assembly) ->
      let ref_db = List.Assoc.find_exn ref_blast_dbs s.ref_species in
      let query = trinity_assembly in
      let query_species= s.species in
      let query_id = s.id in
      let ref_transcriptome = (configuration_dir / ref_transcriptomes s.ref_species) in
      let seq2fam = (configuration_dir / ref_seq_fam_links s.ref_species) in
      let r =
        seq_dispatcher
          ~s2s_tab_by_family:true
          ~query
          ~query_species
          ~query_id
          ~ref_transcriptome
          ~seq2fam
          ~ref_db
          ~threads
      in
      (s, r)
    )

let apytram_orfs_ref_fams_of_apytram_annotated_ref_fams apytram_annotated_ref_fams memory =
  List.map apytram_annotated_ref_fams ~f:(fun (s, f, apytram_result_fasta) ->
      if s.run_transdecoder then
        let pep_min_length = 20 in
        let retain_long_orfs = 150 in
        let filtered_orf = Transdecoder.transdecoder ~descr:("Apytram." ^ f) ~only_top_strand:true ~retain_long_orfs ~pep_min_length ~only_best_orf:true ~threads:1 ~memory apytram_result_fasta in
        (s, f, filtered_orf)
      else
        (s, f, apytram_result_fasta)
    )

let checkfamily ?ref_db ~(input:fasta workflow) ~family ~ref_transcriptome ~seq2fam : fasta workflow =
  let tmp_checkfamily = dest // "tmp" in
  let dest_checkfamily = dest // "sequences.fa" in
  workflow ~version:8 ~descr:("CheckFamily.py:" ^ family ^ " ") [
    mkdir_p tmp_checkfamily;
    cd tmp_checkfamily;
    cmd "CheckFamily.py"  [
      opt "-tmp" ident tmp_checkfamily ;
      opt "-i" dep input ;
      opt "-t" dep ref_transcriptome ;
      opt "-f" string family;
      opt "-t2f" dep seq2fam;
      opt "-o" ident dest_checkfamily;
      option (opt "-d" (fun blast_db -> seq [dep blast_db ; string "/db"])) ref_db;
    ]
  ]
  / selector [ "sequences.fa" ]

let apytram_checked_families_of_orfs_ref_fams apytram_orfs_ref_fams configuration_dir ref_blast_dbs =
  List.map apytram_orfs_ref_fams ~f:(fun (s, f, apytram_orfs_fasta) ->
    let input = apytram_orfs_fasta in
    let ref_transcriptome = configuration_dir / ref_transcriptomes s.ref_species in
    let seq2fam = configuration_dir / ref_seq_fam_links s.ref_species in
    let ref_db = List.Assoc.find_exn ref_blast_dbs s.ref_species in
    let checked_families_fasta = checkfamily ~input ~family:f ~ref_transcriptome ~seq2fam ~ref_db in
    (s, f, checked_families_fasta)
    )

let parse_apytram_results apytram_annotated_ref_fams =
  let config = Bistro.Expr.(
      List.map apytram_annotated_ref_fams ~f:(fun (s, f, w) ->
          seq ~sep:"\t" [ string s.species ; string s.id ; string f ; dep w ]
        )
      |> seq ~sep:"\n"
    )
  in
  workflow ~version:4 ~descr:"Parse_apytram_results.py" ~np:1  [
    cmd "Parse_apytram_results.py" [ file_dump config ; dest ]
  ]


let seq_integrator
    ?realign_ali
    ?resolve_polytomy
    ?(species_to_refine_list = [])
    ~family
    ~trinity_fam_results_dirs
    ~apytram_results_dir
    ~alignment_sp2seq
    alignment
  : _ directory workflow =

  let get_trinity_file_list extension dirs =
    List.map  dirs ~f:(fun (s,dir) ->
        [ dep dir ; string ("/Trinity." ^ s.id ^ "." ^ s.species ^ "." ^ family ^ "." ^ extension) ; string ","]
      )
    |> List.concat
  in

  let get_apytram_file_list extension dir =
    [ dep dir ; string ("/apytram." ^ family ^ "." ^ extension) ; string ","]
  in

  let trinity_fasta_list  =  get_trinity_file_list "fa" trinity_fam_results_dirs in
  let trinity_sp2seq_list  =  get_trinity_file_list "sp2seq.txt" trinity_fam_results_dirs in

  let apytram_fasta  =  get_apytram_file_list "fa" apytram_results_dir in
  let apytram_sp2seq  =  get_apytram_file_list "sp2seq.txt" apytram_results_dir in

  let sp2seq = List.concat [[dep alignment_sp2seq ; string "," ] ; trinity_sp2seq_list ; apytram_sp2seq ]  in
  let fasta = List.concat [trinity_fasta_list; apytram_fasta]  in

  let tmp_merge = dest // "tmp" in

  workflow ~version:10 ~descr:("SeqIntegrator.py:" ^ family ^ " ") [
    mkdir_p tmp_merge ;
    cmd "SeqIntegrator.py"  [
      opt "-tmp" ident tmp_merge;
      opt "-log" seq [ tmp_merge ; string ("/SeqIntegrator." ^ family ^ ".log" )] ;
      opt "-ali" string alignment ;
      opt "-fa" (seq ~sep:"") fasta;
      option (flag string "--realign_ali") realign_ali;
      option (flag string "--resolve_polytomy") resolve_polytomy;
      opt "-sp2seq" (seq ~sep:"") sp2seq  ; (* list de sp2seq delimited by comas *)
      opt "-out" seq [ dest ; string "/" ; string family] ;
      opt "-sptorefine" (seq ~sep:",") (List.map species_to_refine_list ~f:(fun sp -> string sp) );
    ]
  ]


let seq_filter
    ?realign_ali
    ?resolve_polytomy
    ?(species_to_refine_list = [])
    ~filter_threshold
    ~family
    ~alignment
    ~tree
    ~sp2seq
    : _ directory workflow  =

  let tmp_merge = dest // "tmp" in

  workflow ~version:3 ~descr:("SeqFilter.py:" ^ family ^ " ") [
    mkdir_p tmp_merge ;
    cmd "SeqFilter.py"  [
      opt "-tmp" ident tmp_merge;
      opt "-log" seq [ tmp_merge ; string ("/SeqFilter." ^ family ^ ".log" )] ;
      opt "-ali" dep alignment ;
      opt "-t" dep tree;
      opt "--filter_threshold" float filter_threshold;
      option (flag string "--realign_ali") realign_ali;
      option (flag string "--resolve_polytomy") resolve_polytomy;
      opt "-sp2seq" dep sp2seq  ;
      opt "-out" seq [ dest ; string "/" ; string family] ;
      opt "-sptorefine" (seq ~sep:",") (List.map species_to_refine_list ~f:(fun sp -> string sp) );
    ]
  ]


let merged_families_of_families configuration configuration_dir trinity_annotated_fams apytram_results_dir =
  List.map configuration.families ~f:(fun family ->
      let trinity_fam_results_dirs=
        List.map configuration.trinity_samples ~f:(fun s ->
            (s , List.Assoc.find_exn trinity_annotated_fams s)
          )
      in

      let alignment = configuration.alignments_dir ^ "/" ^ family ^ ".fa"  in
      let alignment_sp2seq = configuration_dir / ali_species2seq_links family in
      let species_to_refine_list = List.map configuration.all_ref_samples ~f:(fun s -> s.species) in
      
      let w = seq_integrator ~realign_ali:true ~resolve_polytomy:true ~species_to_refine_list ~family ~trinity_fam_results_dirs ~apytram_results_dir ~alignment_sp2seq  alignment in
      
      let wf = if configuration.ali_sister_threshold > 0. then
                 let filter_threshold = configuration.ali_sister_threshold in
                 let tree = w / selector [family ^ ".tree"] in
                 let alignment = w / selector [family ^ ".fa"] in
                 let sp2seq = w / selector [family ^ ".sp2seq.txt"] in
                 seq_filter ~realign_ali:true ~resolve_polytomy:true ~filter_threshold ~species_to_refine_list ~family ~tree ~alignment ~sp2seq
               else
                 w
               in
      (family, w, wf )
    )




let phyldog_by_fam_of_merged_families merged_families configuration =
  List.map  merged_families ~f:(fun (fam, merged_without_filter_family, merged_and_filtered_family) ->
    let merged_family = if configuration.ali_sister_threshold > 0. then
                            merged_and_filtered_family
                        else
                            merged_without_filter_family
                        in
    let ali = merged_family / selector [ fam ^ ".fa" ] in
    let tree = merged_family / selector [ fam ^ ".tree" ] in
    let link = merged_family / selector [ fam ^ ".sp2seq.txt" ] in
    let sptreefile = configuration.species_tree_file in
    let profileNJ_tree = (ProfileNJ.profileNJ ~sptreefile ~link ~tree ~threshold:1.0 ) / selector [ fam ^ ".tree" ] in
    let threads = Pervasives.min 2 configuration.threads in
    let memory = Pervasives.(configuration.memory / configuration.threads) in
    let topogene = configuration.refinetree in
    (fam, Phyldog.phyldog_by_fam ~threads ~memory ~topogene ~timelimit:9999999 ~sptreefile ~link ~tree:profileNJ_tree ali, merged_family)
    )

let realign_merged_families merged_and_reconciled_families configuration =
  List.map  merged_and_reconciled_families ~f:(fun (fam, reconciled_w, merged_w) ->
    let ali = merged_w / selector [ fam ^ ".fa" ] in
    let treein = reconciled_w / selector [ "Gene_trees/" ^ fam ^ ".ReconciledTree" ] in
    let threads = 1 in
    (fam, Aligner.mafft ~threads ~treein ~auto:false ali, reconciled_w, merged_w)
    )

let merged_families_distributor merged_reconciled_and_realigned_families configuration=
  let extension_list_merged = [(".fa","Merged_fasta");(".tree","Merged_tree");(".sp2seq.txt","Sp2Seq_link")] in
  let extension_list_reconciled = [(".ReconciledTree","Gene_trees/","Reconciled_Gene_tree")] in
  let extension_list_realigned = [(".realign.fa","Realigned_fasta/")] in
  workflow ~descr:"build_output_directory" ~version:1 [
    mkdir_p tmp;

    mkdir_p (dest // "Merged_fasta");
    mkdir_p (dest // "Merged_tree");
    mkdir_p (dest // "Sp2Seq_link");
    if configuration.run_reconciliation then
       mkdir_p (dest // "Reconciled_Gene_tree")
    else
        mkdir_p tmp
    ;
    if configuration.refineali && configuration.run_reconciliation then
      mkdir_p (dest // "Realigned_fasta")
    else
      mkdir_p tmp
    ;
    let script = Bistro.Expr.(
      List.map merged_reconciled_and_realigned_families ~f:(fun (f, realigned_w, reconciled_w, merged_w) ->
          List.concat[
              List.map extension_list_merged ~f:(fun (ext,dir) ->
                let input = merged_w / selector [ f ^ ext ] in
                let output = dest // dir // (f ^ ext)  in
                seq ~sep:" " [ string "ln -s"; dep input ; ident output ]
              )
              ;
              if configuration.run_reconciliation then
                List.concat [
                  List.map extension_list_reconciled ~f:(fun (ext,dirin,dirout) ->
                    let input = reconciled_w / selector [ dirin ^ f ^ ext ] in
                    let output = dest // dirout // (f ^ ext)  in
                    seq ~sep:" " [ string "ln -s"; dep input ; ident output ]
                    )
                  ;
                  if configuration.refineali then
                    List.map extension_list_realigned ~f:(fun (ext,dir) ->
                        let input = realigned_w in
                        let output = dest // dir // (f ^ ext)  in
                        seq ~sep:" " [ string "ln -s"; dep input ; ident output ]
                    )
                  else
                    []
                ;]
              else
                  []
              ;
              ]

            |> seq ~sep:"\n"
          )
        |> seq ~sep:"\n"
      )
    in
    cmd "bash" [ file_dump script ]
  ]

let get_reconstructed_sequences merged_and_reconciled_families_dirs configuration =
  let species_to_refine_list = List.map configuration.all_ref_samples ~f:(fun s -> s.species) in
  workflow ~descr:"GetReconstructedSequences.py" ~version:2 [
    mkdir_p dest;
    cmd "GetReconstructedSequences.py"  [
      dep merged_and_reconciled_families_dirs // "Merged_fasta";
      dep merged_and_reconciled_families_dirs // "Sp2Seq_link";
      seq ~sep:"," (List.map species_to_refine_list ~f:(fun sp -> string sp));
      ident dest
    ]
  ]

let phyldog_of_merged_families_dirs configuration merged_families_dirs =
  let seqdir = merged_families_dirs / selector [ "Merged_fasta" ] in
  let treedir = merged_families_dirs / selector [ "Merged_tree" ] in
  let linkdir = merged_families_dirs / selector [ "Sp2Seq_link" ] in
  let sptreefile = configuration.species_tree_file in
  let threads_max = (List.length configuration.families) + 1 in
  let threads = Pervasives.min threads_max configuration.threads in
  let memory = configuration.memory in
  Phyldog.phyldog ~threads ~memory ~topogene:true ~timelimit:9999999 ~sptreefile ~linkdir ~treedir seqdir



let output_of_phyldog phyldog merged_families families =
  workflow ~descr:"output_of_phyldog" ~version:1 [
    mkdir_p (dest // "Alignments");
    mkdir_p (dest // "Sp2Seq_link");
    mkdir_p (dest // "Gene_trees");
    let extension_list = [(".fa","Alignments");(".sp2seq.txt","Sp2Seq_link")] in
    let script = Bistro.Expr.(
        seq ~sep:"\n" [
          List.map extension_list ~f:(fun (ext,dir) ->
              List.map  merged_families ~f:(fun (f, w) ->
                  let input = w / selector [ f ^ ext ] in
                  let output = dest // dir // (f ^ ext)  in
                  seq ~sep:" " [ string "ln -s"; dep input ; ident output ]
                )
              |> seq ~sep:"\n"
            )
          |> seq ~sep:"\n" ;
          let (ext,dir) = (".ReconciledTree","Gene_trees/") in
          List.map families ~f:(fun f ->
              let input = phyldog / selector [ dir ^ f ^ ext ] in
              let output = dest // dir // (f ^ ".tree")  in
              seq ~sep:" " [ string "ln -s"; dep input ; ident output ]
            )
          |> seq ~sep:"\n";
        ]
      )
    in
    cmd "bash" [ file_dump script ];
  ]



let build_app configuration =

  let divided_memory = Pervasives.(configuration.memory / configuration.threads) in

  let configuration_dir = parse_input configuration in

  let ref_blast_dbs = ref_blast_dbs_of_configuration_dir configuration configuration_dir in

  let fasta_reads = fastq_to_fasta_conversion configuration configuration_dir in

  let norm_fasta = normalize_fasta fasta_reads configuration in

  let trinity_assemblies = trinity_assemblies_of_norm_fasta norm_fasta configuration in

  let trinity_orfs = transdecoder_orfs_of_trinity_assemblies trinity_assemblies configuration in

  let trinity_assemblies_stats = assemblies_stats_of_fasta trinity_assemblies in

  let trinity_orfs_stats = assemblies_stats_of_fasta trinity_orfs in

  let trinity_annotated_fams = trinity_annotated_fams_of_trinity_assemblies configuration_dir ref_blast_dbs configuration trinity_orfs in

  let blast_dbs = blast_dbs_of_norm_fasta norm_fasta in

  let apytram_annotated_ref_fams =
    let pairs = List.cartesian_product configuration.apytram_samples configuration.families in
    List.map pairs ~f:(fun (s, fam) ->
        let query = configuration_dir / ref_fams s.ref_species fam in
        let blast_db = List.Assoc.find_exn blast_dbs s in
        let db_type = sample_fastq_orientation s.sample_fastq in
        let w = Apytram.apytram ~no_best_file:true ~write_even_empty:true ~plot:false ~i:5 ~evalue:1e-5 ~memory:divided_memory ~query db_type blast_db in
        let apytram_filename = "apytram." ^ s.ref_species ^ "." ^ fam ^ ".fasta" in
        (s, fam, w / selector [ apytram_filename ] )
      )
  in

  let apytram_annotated_ref_fams_by_fam =
  List.concat
    (List.map configuration.families ~f:(fun fam ->
    let query = concat (List.map configuration.all_apytram_ref_species ~f:(fun sp -> configuration_dir / ref_fams sp fam)) in
    let blast_dbs = blast_dbs (*(s, w)*) in
    let w = Apytram.apytram_multi_species ~no_best_file:true ~write_even_empty:true ~plot:false ~i:5 ~evalue:1e-5 ~out_by_species:true ~memory:divided_memory ~fam ~query blast_dbs in
    List.map configuration.apytram_samples ~f:(fun s ->
      let apytram_filename = "apytram." ^ fam ^ "." ^ s.id ^ ".fasta" in
      (s, fam, w / selector [ apytram_filename ] )
       )
    )
    )
  in

  let apytram_orfs_ref_fams = apytram_orfs_ref_fams_of_apytram_annotated_ref_fams apytram_annotated_ref_fams_by_fam divided_memory in

  let apytram_checked_families =  apytram_checked_families_of_orfs_ref_fams apytram_orfs_ref_fams configuration_dir ref_blast_dbs in

  let apytram_results_dir = parse_apytram_results apytram_checked_families in

  let merged_families = merged_families_of_families configuration configuration_dir trinity_annotated_fams apytram_results_dir in

  let merged_and_reconciled_families = phyldog_by_fam_of_merged_families merged_families configuration in

  let merged_reconciled_and_realigned_families = realign_merged_families merged_and_reconciled_families configuration in

  let merged_reconciled_and_realigned_families_dirs = merged_families_distributor merged_reconciled_and_realigned_families configuration in

  let reconstructed_sequences = get_reconstructed_sequences merged_reconciled_and_realigned_families_dirs configuration in

  (*let phyldog = phyldog_of_merged_families_dirs configuration merged_families_dirs in

  let output = output_of_phyldog phyldog merged_families configuration.families in
*)


  let open Bistro_app in

  let target_to_sample_fasta s d = function
    | Fasta_Single_end (w, _ ) -> [[ d ; s.id ^ "_" ^ s.species ^ ".fa" ] %> w ]
    | Fasta_Paired_end (lw, rw , _) -> [[ d ; s.id ^ "_" ^ s.species ^ ".left.fa" ] %> lw ; [ d ; s.id ^ "_" ^ s.species ^ ".right.fa" ] %> lw]
  in
  let repo =
    List.concat [
      if configuration.debug then
      List.concat [
        [[ "tmp" ; "configuration" ] %>  configuration_dir ]
        ;
        List.concat (List.map fasta_reads ~f:(fun (s,sample_fasta) -> target_to_sample_fasta s "tmp/rna_seq/raw_fasta" sample_fasta))
        ;
        List.concat (List.map norm_fasta ~f:(fun (s,norm_fasta) -> target_to_sample_fasta s "tmp/rna_seq/norm_fasta" norm_fasta))
        ;
        List.map trinity_assemblies ~f:(fun (s,trinity_assembly) ->
            [ "tmp" ; "trinity_assembly" ; "trinity_assemblies" ; "Trinity_assemblies." ^ s.id ^ "_" ^ s.species ^ ".fa" ] %> trinity_assembly
          )
        ;
        List.map trinity_orfs ~f:(fun (s,trinity_orf) ->
            [ "tmp" ; "trinity_assembly" ; "trinity_assemblies" ; "Transdecoder_cds." ^ s.id ^ "_" ^ s.species ^ ".fa" ] %> trinity_orf
          )
        ;
        List.map trinity_assemblies_stats ~f:(fun (s,trinity_assembly_stats) ->
            [ "tmp" ; "trinity_assembly" ; "trinity_assemblies_stats" ; "Trinity_assemblies." ^ s.id ^ "_" ^ s.species ^ ".stats" ] %> trinity_assembly_stats
          )
        ;
        List.map trinity_orfs_stats ~f:(fun (s,trinity_orfs_stats) ->
            [ "tmp" ; "trinity_assembly" ; "trinity_assemblies_stats" ; "Transdecoder_cds." ^ s.id ^ "_" ^ s.species ^ ".stats" ] %> trinity_orfs_stats
          )
        ;
        List.map trinity_annotated_fams ~f:(fun (s,trinity_annotated_fams) ->
            [ "tmp" ; "trinity_blast_annotation" ; "trinity_annotated_fams" ; s.id ^ "_" ^ s.species ^ ".vs." ^ s.ref_species ] %> trinity_annotated_fams
          )
        ;
         List.map ref_blast_dbs ~f:(fun (ref_species, blast_db) ->
            [ "tmp" ; "trinity_blast_annotation" ; "ref_blast_db" ; ref_species ] %> blast_db
          )
        ;
        List.map blast_dbs ~f:(fun (s,blast_db) ->
            [ "tmp" ; "rna_seq" ;"blast_db" ; s.id ^ "_" ^ s.species ] %> blast_db
          )
        ;
       (* List.map apytram_annotated_ref_fams ~f:(fun (s, fam, apytram_result) ->
            [ "tmp" ; "apytram_assembly" ; "apytram_annotated_fams" ; fam ; s.id ^ "_" ^ s.species ^ ".fa" ] %> apytram_result
          )
        ;
        *)
        List.map apytram_annotated_ref_fams_by_fam ~f:(fun (s, fam, apytram_result) ->
            [ "tmp" ; "apytram_assembly" ; "apytram_annotated_fams_by_fam" ; fam ; s.id ^ "_" ^ s.species ^ ".fa" ] %> apytram_result
          )
        ;
        List.map apytram_orfs_ref_fams ~f:(fun (s, fam, apytram_result) ->
            [ "tmp" ; "apytram_assembly" ; "apytram_transdecoder_orfs" ; fam ; s.id ^ "_" ^ s.species ^ ".fa" ] %> apytram_result
          )
        ;
        List.map apytram_checked_families ~f:(fun (s, fam, apytram_result) ->
            [ "tmp" ; "apytram_assembly" ; "apytram_checked_families" ; fam ; s.id ^ "_" ^ s.species ^ ".fa"] %> apytram_result
          )
        ;
        [["tmp" ; "apytram_assembly" ;"apytram_results" ] %> apytram_results_dir]
        ;
        List.map merged_families ~f:(fun (fam, merged_family, merged_and_filtered_family) ->
            [ "tmp" ; "merged_families" ; fam  ] %> merged_family;
            
          )
        ;
        List.map merged_families ~f:(fun (fam, merged_family, merged_and_filtered_family) ->
            [ "tmp" ; "merged_filtered_families" ; fam  ] %> merged_and_filtered_family;
            
          )
      ]
      else
      []
      ;
      [["merged_families_dir"] %> merged_reconciled_and_realigned_families_dirs]
      ;
      [["reconstructed_sequences"] %> reconstructed_sequences]
      ;
    ]
  in
  let repo_app = Bistro_app.of_repo repo ~outdir:configuration.outdir in
  let stats_app =
    List.map trinity_assemblies_stats ~f:(fun (s, trinity_assembly_stats) -> (s, pureW trinity_assembly_stats))
    |> assoc
  in
  let f_app = pure (fun () trinity_assemblies_stats -> Report.generate ~trinity_assemblies_stats (Filename.concat configuration.outdir "report.html")) in
  f_app $ repo_app $ stats_app
