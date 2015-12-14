import os
import sys
import logging
import subprocess

## TO DO: check if Trinity is installed and is path is Trinity v2.2>=

class Trinity:
    """Define an object to launch Trinity"""
    def __init__(self, InputFile, OutputFile):
        self.seqType = "fa"
        self.max_memory = 1
        self.CPU = 1
        self.InputFile = InputFile
        self.OutputFile = OutputFile
        self.MinLength = 200
        self.Paired = False
        self.FullCleanup = False
        self.NoPathMerging = False

    def launch(self):
        ExitCode = 0
        command = ["Trinity","--no_version_check","--seqType", self.seqType,"--single", self.InputFile,
                   "--output", self.OutputFile,
                   "--CPU", str(self.CPU), "--max_memory", str(self.max_memory)+"G"]
       
        if self.FullCleanup:
            command.append("--full_cleanup")
            
        if self.MinLength != 200:
            command.extend(["--min_contig_length",str(self.MinLength)])
        
        if self.Paired:
            command.append("--run_as_paired")
        
        if self.NoPathMerging:
            command.append("--no_path_merging")
            
        try:
            ExitCode = subprocess.call(command,
                                       stdout=open("/dev/null", "w"),
                                       stderr=open("/dev/null", "w"))
        except:
            os.system("echo Unexpected error when we launch Trinity:\n")
            print " ".join(command)
            
        return ExitCode
    
###############################################################################
#
#     ______  ____   ____  ____   ____  ______  __ __
#    |      ||    \ |    ||    \ |    ||      ||  |  |
#    |      ||  D  ) |  | |  _  | |  | |      ||  |  |
#    |_|  |_||    /  |  | |  |  | |  | |_|  |_||  ~  |
#      |  |  |    \  |  | |  |  | |  |   |  |  |___, |
#      |  |  |  .  \ |  | |  |  | |  |   |  |  |     |
#      |__|  |__|\_||____||__|__||____|  |__|  |____/
#
###############################################################################
#
# Required:
#
#  --seqType <string>      :type of reads: ( fa, or fq )
#
#  --max_memory <string>      :suggested max memory to use by Trinity where limiting can be enabled. (jellyfish, sorting, etc)
#                            provied in Gb of RAM, ie.  '--max_memory 10G'
#
#  If paired reads:
#      --left  <string>    :left reads, one or more file names (separated by commas, no spaces)
#      --right <string>    :right reads, one or more file names (separated by commas, no spaces)
#
#  Or, if unpaired reads:
#      --single <string>   :single reads, one or more file names, comma-delimited (note, if single file contains pairs, can use flag: --run_as_paired )
#
####################################
##  Misc:  #########################
#
#  --SS_lib_type <string>          :Strand-specific RNA-Seq read orientation.
#                                   if paired: RF or FR,
#                                   if single: F or R.   (dUTP method = RF)
#                                   See web documentation.
#
#  --CPU <int>                     :number of CPUs to use, default: 2
#  --min_contig_length <int>       :minimum assembled contig length to report
#                                   (def=200)
#
#  --long_reads <string>           :fasta file containing error-corrected or circular consensus (CCS) pac bio reads
#
#  --genome_guided_bam <string>    :genome guided mode, provide path to coordinate-sorted bam file.
#                                   (see genome-guided param section under --show_full_usage_info)
#
#  --jaccard_clip                  :option, set if you have paired reads and
#                                   you expect high gene density with UTR
#                                   overlap (use FASTQ input file format
#                                   for reads).
#                                   (note: jaccard_clip is an expensive
#                                   operation, so avoid using it unless
#                                   necessary due to finding excessive fusion
#                                   transcripts w/o it.)
#
#  --trimmomatic                   :run Trimmomatic to quality trim reads
#                                        see '--quality_trimming_params' under full usage info for tailored settings.
#                                  
#
#  --normalize_reads               :run in silico normalization of reads. Defaults to max. read coverage of 50.
#                                       see '--normalize_max_read_cov' under full usage info for tailored settings.
#     
#  --no_distributed_trinity_exec   :do not run Trinity phase 2 (assembly of partitioned reads), and stop after generating command list.
#
#
#  --output <string>               :name of directory for output (will be
#                                   created if it doesn't already exist)
#                                   default( your current working directory: "/home/crey02/Documents/Projets/Convergences/Pipeline/apytram/trinity_out_dir" 
#                                    note: must include 'trinity' in the name as a safety precaution! )
#  
#  --full_cleanup                  :only retain the Trinity fasta file, rename as ${output_dir}.Trinity.fasta
#
#  --cite                          :show the Trinity literature citation
#
#  --verbose                       :provide additional job status info during the run.
#
#  --version                       :reports Trinity version (v2.1.1) and exits.
#
#  --show_full_usage_info          :show the many many more options available for running Trinity (expert usage).

#
#  --KMER_SIZE <int>               :kmer length to use (default: 25)  max=32
#
#  --prep                          :Only prepare files (high I/O usage) and stop before kmer counting.
#
#  --no_cleanup                    :retain all intermediate input files.
#
#  --no_version_check              :dont run a network check to determine if software updates are available.
#
####################################################
# Inchworm and K-mer counting-related options: #####
#
#  --min_kmer_cov <int>           :min count for K-mers to be assembled by
#                                  Inchworm (default: 1)
#  --inchworm_cpu <int>           :number of CPUs to use for Inchworm, default is min(6, --CPU option)
#
#  --no_run_inchworm              :stop after running jellyfish, before inchworm. (phase 1, read clustering only)
#
###################################
# Chrysalis-related options: ######
#
#  --max_reads_per_graph <int>    :maximum number of reads to anchor within
#                                  a single graph (default: 200000)
#  --min_glue <int>               :min number of reads needed to glue two inchworm contigs
#                                  together. (default: 2) 
#
#  --no_bowtie                    :dont run bowtie to use pair info in chrysalis clustering.
#
#  --no_run_chrysalis             :stop after running inchworm, before chrysalis. (phase 1, read clustering only)
#
#####################################
###  Butterfly-related options:  ####
#
#  --bfly_opts <string>            :additional parameters to pass through to butterfly
#                                   (see butterfly options: java -jar Butterfly.jar ).
#                                   (note: only for expert or experimental use.  Commonly used parameters are exposed through this Trinity menu here).
#
#    //////////////////////////////////
#    Alternative reconstruction modes:
#                                  Default mode is the 'regular' Butterfly transcript reconstruction by graph node extension.
#
#       --PasaFly                  PASA-like algorithm for maximally-supported isoforms 
#           or
#       --CuffFly                  Cufflinks-like algorithm to report minimum transcripts
#
#
#  Butterfly read-pair grouping settings (used for all reconstruction modes to define 'pair paths'):
#
#  --group_pairs_distance <int>    :maximum length expected between fragment pairs (default: 500)
#                                   (reads outside this distance are treated as single-end)
#
#  ///////////////////////////////////////////////
#  Butterfly default reconstruction mode settings. (no CuffFly or PasaFly custom settings are currently available).
#                                   
#  --path_reinforcement_distance <int>   :minimum overlap of reads with growing transcript 
#                                         path (default: PE: 75, SE: 25)
#                                         Set to 1 for the most lenient path extension requirements.
#
#
#  /////////////////////////////////////////
#  Butterfly transcript reduction settings:
#
#  --no_path_merging            : all final transcript candidates are output (including SNP variations, however, some SNPs may be unphased)  
#
#  By default, alternative transcript candidates are merged (in reality, discarded) if they are found to be too similar, according to the following logic:
#
#  (identity=(numberOfMatches/shorterLen) > 95.0% or if we have <= 2 mismatches) and if we have internal gap lengths <= 10
#
#  with parameters as:
#      
#      --min_per_id_same_path <int>          default: 98     min percent identity for two paths to be merged into single paths
#      --max_diffs_same_path <int>           default: 2      max allowed differences encountered between path sequences to combine them
#      --max_internal_gap_same_path <int>    default: 10     maximum number of internal consecutive gap characters allowed for paths to be merged into single paths.
#
#      If, in a comparison between two alternative transcripts, they are found too similar, the transcript with the greatest cumulative 
#      compatible read (pair-path) support is retained, and the other is discarded.
#
#
#  //////////////////////////////////////////////
#  Butterfly Java and parallel execution settings.
#
#  --bflyHeapSpaceMax <string>     :java max heap space setting for butterfly
#                                   (default: 4G) => yields command
#                  'java -Xmx4G -jar Butterfly.jar ... $bfly_opts'
#  --bflyHeapSpaceInit <string>    :java initial hap space settings for
#                                   butterfly (default: 1G) => yields command
#                  'java -Xms1G -jar Butterfly.jar ... $bfly_opts'
#  --bflyGCThreads <int>           :threads for garbage collection
#                                   (default: 2))
#  --bflyCPU <int>                 :CPUs to use (default will be normal 
#                                   number of CPUs; e.g., 2)
#  --bflyCalculateCPU              :Calculate CPUs based on 80% of max_memory
#                                   divided by maxbflyHeapSpaceMax
#
#  --bfly_jar <string>             : /path/to/Butterfly.jar, otherwise default
#                                    Trinity-installed version is used. 
#                                    
#
#
################################################################################
#### Quality Trimming Options ####  
# 
#  --quality_trimming_params <string>   defaults to: "ILLUMINACLIP:/usr/local/bin/trinityrnaseq-2.1.1/trinity-plugins/Trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25"
#
################################################################################
####  In silico Read Normalization Options ###
#
#  --normalize_max_read_cov <int>       defaults to 50
#  --normalize_by_read_set              run normalization separate for each pair of fastq files,
#                                       then one final normalization that combines the individual normalized reads.
#                                       Consider using this if RAM limitations are a consideration.
#
################################################################################
#### Genome-guided de novo assembly
# 
#  * required:
#
# --genome_guided_max_intron <int>     :maximum allowed intron length (also maximum fragment span on genome)
#
#  * optional:
#
# --genome_guided_min_coverage <int>   :minimum read coverage for identifying and expressed region of the genome. (default: 1)
#
# --genome_guided_min_reads_per_partition <int>   :default min of 10 reads per partition
#
#
#################################
# Grid-computing options: #######
#
#  --grid_conf <string>                 :configuration file for supported compute farms
#                                       ex.  TRINITY_HOME/htc_conf/BroadInst_LSF.conf
#                                       currently supported computing gris: LSF, SGE
#
#  --grid_node_CPU <int>                number of threads for each parallel process to leverage. (default: 1)
#
#  --grid_node_max_memory <string>         max memory targeted for each grid node. (default: 1G)
#
#
    #
#
###############################################################################
#
#  *Note, a typical Trinity command might be:
#
#        Trinity --seqType fq --max_memory 50G --left reads_1.fq  --right reads_2.fq --CPU 6
#
#
#    and for Genome-guided Trinity:
#
#        Trinity --genome_guided_bam rnaseq_alignments.csorted.bam --max_memory 50G
#                --genome_guided_max_intron 10000 --CPU 6
#
#     see: /usr/local/bin/trinityrnaseq-2.1.1/sample_data/test_Trinity_Assembly/
#          for sample data and 'runMe.sh' for example Trinity execution
#
#     For more details, visit: http://trinityrnaseq.github.io
#
###############################################################################
