import os
import sys
import logging
import subprocess


class Makeblastdb:
    """Define an object to create a local database"""
    def __init__(self, InputFile, OutputFiles):
        self.logger = logging.getLogger('apytram.lib.BlastPlus.Makeblastdb')
        self.Title = ""
        self.InputFile = InputFile
        self.OutputFiles = OutputFiles
        self.LogFile = ""
        self.Dbtype = "nucl"
        self.IndexedDatabase = True 

    def launch(self):
        ExitCode = 1
        command = ["makeblastdb","-in",self.InputFile,"-out", os.path.abspath(self.OutputFiles),
                   "-dbtype", self.Dbtype ]
        self.logger.debug(" ".join(command))
        if self.IndexedDatabase:
            command.append("-parse_seqids")
        try:
            ExitCode = subprocess.call(command,
                                       stdout=open("/dev/null", "w"),
                                       stderr=open("/dev/null", "w"))
        except:
            os.system("echo Unexpected error: "+" ".join(command)+"\n")
        return ExitCode
    
#  makeblastdb [-h] [-help] [-in input_file] [-input_type type]
#    -dbtype molecule_type [-title database_title] [-parse_seqids]
#    [-hash_index] [-mask_data mask_data_files] [-gi_mask]
#    [-gi_mask_name gi_based_mask_names] [-out database_name]
#    [-max_file_sz number_of_bytes] [-taxid TaxID] [-taxid_map TaxIDMapFile]
#    [-logfile File_Name] [-version]
#
#DESCRIPTION
#   Application to create BLAST databases, version 2.2.28+
#
#REQUIRED ARGUMENTS
# -dbtype <String, `nucl', `prot'>
#   Molecule type of target db
#
#OPTIONAL ARGUMENTS
# -h
#   Print USAGE and DESCRIPTION;  ignore all other parameters
# -help
#   Print USAGE, DESCRIPTION and ARGUMENTS; ignore all other parameters
# -version
#   Print version number;  ignore other arguments
#
# *** Input options
# -in <File_In>
#   Input file/database name
#   Default = `-'
# -input_type <String, `asn1_bin', `asn1_txt', `blastdb', `fasta'>
#   Type of the data specified in input_file
#   Default = `fasta'
#
# *** Configuration options
# -title <String>
#   Title for BLAST database
#   Default = input file name provided to -in argument
# -parse_seqids
#   Option to parse seqid for FASTA input if set, for all other input types
#   seqids are parsed automatically
# -hash_index
#   Create index of sequence hash values.
#
# *** Sequence masking options
# -mask_data <String>
#   Comma-separated list of input files containing masking data as produced by
#   NCBI masking applications (e.g. dustmasker, segmasker, windowmasker)
# -gi_mask
#   Create GI indexed masking data.
#    * Requires:  parse_seqids
# -gi_mask_name <String>
#   Comma-separated list of masking data output files.
#    * Requires:  mask_data, gi_mask
#
# *** Output options
# -out <String>
#   Name of BLAST database to be created
#   Default = input file name provided to -in argumentRequired if multiple
#   file(s)/database(s) are provided as input
# -max_file_sz <String>
#   Maximum file size for BLAST database files
#   Default = `1GB'
#
# *** Taxonomy options
# -taxid <Integer, >=0>
#   Taxonomy ID to assign to all sequences
#    * Incompatible with:  taxid_map
# -taxid_map <File_In>
#   Text file mapping sequence IDs to taxonomy IDs.
#   Format:<SequenceId> <TaxonomyId><newline>
#    * Incompatible with:  taxid
# -logfile <File_Out>
#   File to which the program log should be redirected


class Blast:
    """Define a object to lauch a blast on a local database"""
    def __init__(self, Program, Database, QueryFile):
        self.logger = logging.getLogger('main.lib.BlastPlus.Blast')
        self.Program = Program
        self.Database = Database
        self.QueryFile = QueryFile
        self.OutputFile = ""
        self.Threads = 1
        self.OutFormat = 8
        self.Evalue = 10
        self.max_target_seqs = 1000000000
        self.perc_identity = 0
        self.Task = ""
   
    def launch(self,OutputFile):
        ExitCode = 1
        if self.Program in ["blastn","blastx","tblastn","tblastx"]:
            command = [self.Program,"-db", self.Database, "-query" , self.QueryFile,
                      "-evalue", str(self.Evalue), "-outfmt", str(self.OutFormat),
                      "-out", OutputFile,
                      "-max_target_seqs", str(self.max_target_seqs),
                      "-num_threads", str(self.Threads)]
            
            if self.perc_identity:
                command.extend(["-perc_identity" , str(self.perc_identity)])
            if self.Task:
                command.extend(["-task",self.Task])
                
            self.logger.debug(" ".join(command))
            try:
                 ExitCode = subprocess.call(command,
                                       stdout=open("/dev/null", "w"),
                                       stderr=open("/dev/null", "w"))
            except:
                os.system("echo Unexpected error\n")
                print " ".join(command)
            return ExitCode
        else:
                self.logger.error("%s not in [blastn,blastx,tblastn,tblastx]" % self.Program)
                return ExitCode
    
#USAGE
#  blastn [-h] [-help] [-import_search_strategy filename]
#    [-export_search_strategy filename] [-task task_name] [-db database_name]
#    [-dbsize num_letters] [-gilist filename] [-seqidlist filename]
#    [-negative_gilist filename] [-entrez_query entrez_query]
#    [-db_soft_mask filtering_algorithm] [-db_hard_mask filtering_algorithm]
#    [-subject subject_input_file] [-subject_loc range] [-query input_file]
#    [-out output_file] [-evalue evalue] [-word_size int_value]
#    [-gapopen open_penalty] [-gapextend extend_penalty]
#    [-perc_identity float_value] [-xdrop_ungap float_value]
#    [-xdrop_gap float_value] [-xdrop_gap_final float_value]
#    [-searchsp int_value] [-max_hsps_per_subject int_value] [-penalty penalty]
#    [-reward reward] [-no_greedy] [-min_raw_gapped_score int_value]
#    [-template_type type] [-template_length int_value] [-dust DUST_options]
#    [-filtering_db filtering_database]
#    [-window_masker_taxid window_masker_taxid]
#    [-window_masker_db window_masker_db] [-soft_masking soft_masking]
#    [-ungapped] [-culling_limit int_value] [-best_hit_overhang float_value]
#    [-best_hit_score_edge float_value] [-window_size int_value]
#    [-off_diagonal_range int_value] [-use_index boolean] [-index_name string]
#    [-lcase_masking] [-query_loc range] [-strand strand] [-parse_deflines]
#    [-outfmt format] [-show_gis] [-num_descriptions int_value]
#    [-num_alignments int_value] [-html] [-max_target_seqs num_sequences]
#    [-num_threads int_value] [-remote] [-version]
#
#DESCRIPTION
#   Nucleotide-Nucleotide BLAST 2.2.28+
#
#OPTIONAL ARGUMENTS
# -h
#   Print USAGE and DESCRIPTION;  ignore all other parameters
# -help
#   Print USAGE, DESCRIPTION and ARGUMENTS; ignore all other parameters
# -version
#   Print version number;  ignore other arguments
#
# *** Input query options
# -query <File_In>
#   Input file name
#   Default = `-'
# -query_loc <String>
#   Location on the query sequence in 1-based offsets (Format: start-stop)
# -strand <String, `both', `minus', `plus'>
#   Query strand(s) to search against database/subject
#   Default = `both'
#
# *** General search options
# -task <String, Permissible values: 'blastn' 'blastn-short' 'dc-megablast'
#                'megablast' 'rmblastn' >
#   Task to execute
#   Default = `megablast'
# -db <String>
#   BLAST database name
#    * Incompatible with:  subject, subject_loc
# -out <File_Out>
#   Output file name
#   Default = `-'
# -evalue <Real>
#   Expectation value (E) threshold for saving hits 
#   Default = `10'
# -word_size <Integer, >=4>
#   Word size for wordfinder algorithm (length of best perfect match)
# -gapopen <Integer>
#   Cost to open a gap
# -gapextend <Integer>
#   Cost to extend a gap
# -penalty <Integer, <=0>
#   Penalty for a nucleotide mismatch
# -reward <Integer, >=0>
#   Reward for a nucleotide match
# -use_index <Boolean>
#   Use MegaBLAST database index
# -index_name <String>
#   MegaBLAST database index name
#
# *** BLAST-2-Sequences options
# -subject <File_In>
#   Subject sequence(s) to search
#    * Incompatible with:  db, gilist, seqidlist, negative_gilist,
#   db_soft_mask, db_hard_mask
# -subject_loc <String>
#   Location on the subject sequence in 1-based offsets (Format: start-stop)
#    * Incompatible with:  db, gilist, seqidlist, negative_gilist,
#   db_soft_mask, db_hard_mask, remote
#
# *** Formatting options
# -outfmt <String>
#   alignment view options:
#     0 = pairwise,
#     1 = query-anchored showing identities,
#     2 = query-anchored no identities,
#     3 = flat query-anchored, show identities,
#     4 = flat query-anchored, no identities,
#     5 = XML Blast output,
#     6 = tabular,
#     7 = tabular with comment lines,
#     8 = Text ASN.1,
#     9 = Binary ASN.1,
#    10 = Comma-separated values,
#    11 = BLAST archive format (ASN.1) 
#   
#   Options 6, 7, and 10 can be additionally configured to produce
#   a custom format specified by space delimited format specifiers.
#   The supported format specifiers are:
#   	    qseqid means Query Seq-id
#   	       qgi means Query GI
#   	      qacc means Query accesion
#   	   qaccver means Query accesion.version
#   	      qlen means Query sequence length
#   	    sseqid means Subject Seq-id
#   	 sallseqid means All subject Seq-id(s), separated by a ';'
#   	       sgi means Subject GI
#   	    sallgi means All subject GIs
#   	      sacc means Subject accession
#   	   saccver means Subject accession.version
#   	   sallacc means All subject accessions
#   	      slen means Subject sequence length
#   	    qstart means Start of alignment in query
#   	      qend means End of alignment in query
#   	    sstart means Start of alignment in subject
#   	      send means End of alignment in subject
#   	      qseq means Aligned part of query sequence
#   	      sseq means Aligned part of subject sequence
#   	    evalue means Expect value
#   	  bitscore means Bit score
#   	     score means Raw score
#   	    length means Alignment length
#   	    pident means Percentage of identical matches
#   	    nident means Number of identical matches
#   	  mismatch means Number of mismatches
#   	  positive means Number of positive-scoring matches
#   	   gapopen means Number of gap openings
#   	      gaps means Total number of gaps
#   	      ppos means Percentage of positive-scoring matches
#   	    frames means Query and subject frames separated by a '/'
#   	    qframe means Query frame
#   	    sframe means Subject frame
#   	      btop means Blast traceback operations (BTOP)
#   	   staxids means Subject Taxonomy ID(s), separated by a ';'
#   	 sscinames means Subject Scientific Name(s), separated by a ';'
#   	 scomnames means Subject Common Name(s), separated by a ';'
#   	sblastnames means Subject Blast Name(s), separated by a ';'
#   			 (in alphabetical order)
#   	sskingdoms means Subject Super Kingdom(s), separated by a ';'
#   			 (in alphabetical order) 
#   	    stitle means Subject Title
#   	salltitles means All Subject Title(s), separated by a '<>'
#   	   sstrand means Subject Strand
#   	     qcovs means Query Coverage Per Subject
#   	   qcovhsp means Query Coverage Per HSP
#   When not provided, the default value is:
#   'qseqid sseqid pident length mismatch gapopen qstart qend sstart send
#   evalue bitscore', which is equivalent to the keyword 'std'
#   Default = `0'
# -show_gis
#   Show NCBI GIs in deflines?
# -num_descriptions <Integer, >=0>
#   Number of database sequences to show one-line descriptions for
#   Not applicable for outfmt > 4
#   Default = `500'
#    * Incompatible with:  max_target_seqs
# -num_alignments <Integer, >=0>
#   Number of database sequences to show alignments for
#   Default = `250'
#    * Incompatible with:  max_target_seqs
# -html
#   Produce HTML output?
#
# *** Query filtering options
# -dust <String>
#   Filter query sequence with DUST (Format: 'yes', 'level window linker', or
#   'no' to disable)
#   Default = `20 64 1'
# -filtering_db <String>
#   BLAST database containing filtering elements (i.e.: repeats)
# -window_masker_taxid <Integer>
#   Enable WindowMasker filtering using a Taxonomic ID
# -window_masker_db <String>
#   Enable WindowMasker filtering using this repeats database.
# -soft_masking <Boolean>
#   Apply filtering locations as soft masks
#   Default = `true'
# -lcase_masking
#   Use lower case filtering in query and subject sequence(s)?
#
# *** Restrict search or results
# -gilist <String>
#   Restrict search of database to list of GI's
#    * Incompatible with:  negative_gilist, seqidlist, remote, subject,
#   subject_loc
# -seqidlist <String>
#   Restrict search of database to list of SeqId's
#    * Incompatible with:  gilist, negative_gilist, remote, subject,
#   subject_loc
# -negative_gilist <String>
#   Restrict search of database to everything except the listed GIs
#    * Incompatible with:  gilist, seqidlist, remote, subject, subject_loc
# -entrez_query <String>
#   Restrict search with the given Entrez query
#    * Requires:  remote
# -db_soft_mask <String>
#   Filtering algorithm ID to apply to the BLAST database as soft masking
#    * Incompatible with:  db_hard_mask, subject, subject_loc
# -db_hard_mask <String>
#   Filtering algorithm ID to apply to the BLAST database as hard masking
#    * Incompatible with:  db_soft_mask, subject, subject_loc
# -perc_identity <Real, 0..100>
#   Percent identity
# -culling_limit <Integer, >=0>
#   If the query range of a hit is enveloped by that of at least this many
#   higher-scoring hits, delete the hit
#    * Incompatible with:  best_hit_overhang, best_hit_score_edge
# -best_hit_overhang <Real, (>=0 and =<0.5)>
#   Best Hit algorithm overhang value (recommended value: 0.1)
#    * Incompatible with:  culling_limit
# -best_hit_score_edge <Real, (>=0 and =<0.5)>
#   Best Hit algorithm score edge value (recommended value: 0.1)
#    * Incompatible with:  culling_limit
# -max_target_seqs <Integer, >=1>
#   Maximum number of aligned sequences to keep 
#   Not applicable for outfmt <= 4
#   Default = `500'
#    * Incompatible with:  num_descriptions, num_alignments
#
# *** Discontiguous MegaBLAST options
# -template_type <String, `coding', `coding_and_optimal', `optimal'>
#   Discontiguous MegaBLAST template type
#    * Requires:  template_length
# -template_length <Integer, Permissible values: '16' '18' '21' >
#   Discontiguous MegaBLAST template length
#    * Requires:  template_type
#
# *** Statistical options
# -dbsize <Int8>
#   Effective length of the database 
# -searchsp <Int8, >=0>
#   Effective length of the search space
# -max_hsps_per_subject <Integer, >=0>
#   Override maximum number of HSPs per subject to save for ungapped searches
#   (0 means do not override)
#   Default = `0'
#
# *** Search strategy options
# -import_search_strategy <File_In>
#   Search strategy to use
#    * Incompatible with:  export_search_strategy
# -export_search_strategy <File_Out>
#   File name to record the search strategy used
#    * Incompatible with:  import_search_strategy
#
# *** Extension options
# -xdrop_ungap <Real>
#   X-dropoff value (in bits) for ungapped extensions
# -xdrop_gap <Real>
#   X-dropoff value (in bits) for preliminary gapped extensions
# -xdrop_gap_final <Real>
#   X-dropoff value (in bits) for final gapped alignment
# -no_greedy
#   Use non-greedy dynamic programming extension
# -min_raw_gapped_score <Integer>
#   Minimum raw gapped score to keep an alignment in the preliminary gapped and
#   traceback stages
# -ungapped
#   Perform ungapped alignment only?
# -window_size <Integer, >=0>
#   Multiple hits window size, use 0 to specify 1-hit algorithm
# -off_diagonal_range <Integer, >=0>
#   Number of off-diagonals to search for the 2nd hit, use 0 to turn off
#   Default = `0'
#
# *** Miscellaneous options
# -parse_deflines
#   Should the query and subject defline(s) be parsed?
# -num_threads <Integer, >=1>
#   Number of threads (CPUs) to use in the BLAST search
#   Default = `1'
#    * Incompatible with:  remote
# -remote
#   Execute search remotely?
#    * Incompatible with:  gilist, seqidlist, negative_gilist, subject_loc,
#   num_threads
#
class Blastdbcmd:
    """Define a object to launch blastdbcmd on a local database"""
    def __init__(self, Database, SequenceNamesFile):
        self.logger = logging.getLogger('main.lib.BlastPlus.Blastdbcmd')
        self.InputFile = SequenceNamesFile
        self.Database = Database
        self.OutputFile = ""
        self.Dbtype = "nucl"

    def launch(self):
		command = ["blastdbcmd", "-db", self.Database, "-entry_batch", self.InputFile, "-dbtype", self.Dbtype , "-out", self.OutputFile]
		self.logger.debug(" ".join(command))
		try:
			BashProcess = subprocess.Popen(command)
		except:
			self.logger.error("Unexpected error when we launched blastdbcmd\n")
		return command
		
    def get_output(self, output = ""):
        Out = ""
        command = ["blastdbcmd", "-db", self.Database, "-entry_batch", self.InputFile, "-dbtype", self.Dbtype ]
		
        if output or self.OutputFile:
            if output:
                self.OutputFile = output
            command.extend(["-out",self.OutputFile])
        
        self.logger.debug(" ".join(command))
        try:
            Out = subprocess.check_output(command,
                                          stderr=open("/dev/null", "w"))
        except:
            self.error.debug("Unexpected error when we launch Blastdbcmd:\n")
            
        return Out

    
    def is_database(self):
        Out = False
        command = ["blastdbcmd","-db",self.Database,"-info"]
        ExitCode = subprocess.call(command,
                                   stdout=open("/dev/null", "w"),
                                   stderr=open("/dev/null", "w"))
        # If no error the database exist
        if not ExitCode:
            Out = True             
        return Out
        
        
#USAGE
#  blastdbcmd [-h] [-help] [-db dbname] [-dbtype molecule_type]
#    [-entry sequence_identifier] [-entry_batch input_file] [-pig PIG] [-info]
#    [-range numbers] [-strand strand] [-mask_sequence_with mask_algo_id]
#    [-out output_file] [-outfmt format] [-target_only] [-get_dups]
#    [-line_length number] [-ctrl_a] [-show_blastdb_search_path]
#    [-list directory] [-remove_redundant_dbs] [-recursive]
#    [-list_outfmt format] [-logfile File_Name] [-version]
#
#DESCRIPTION
#   BLAST database client, version 2.2.28+
#
#OPTIONAL ARGUMENTS
# -h
#   Print USAGE and DESCRIPTION;  ignore all other parameters
# -help
#   Print USAGE, DESCRIPTION and ARGUMENTS; ignore all other parameters
# -version
#   Print version number;  ignore other arguments
#
# *** BLAST database options
# -db <String>
#   BLAST database name
#   Default = `nr'
#    * Incompatible with:  list, recursive, remove_redundant_dbs, list_outfmt,
#   show_blastdb_search_path
# -dbtype <String, `guess', `nucl', `prot'>
#   Molecule type stored in BLAST database
#   Default = `guess'
#
# *** Retrieval options
# -entry <String>
#   Comma-delimited search string(s) of sequence identifiers:
#   	e.g.: 555, AC147927, 'gnl|dbname|tag', or 'all' to select all
#   	sequences in the database
#    * Incompatible with:  entry_batch, pig, info, list, recursive,
#   remove_redundant_dbs, list_outfmt, show_blastdb_search_path
# -entry_batch <File_In>
#   Input file for batch processing (Format: one entry per line, seq id 
#   followed by optional space-delimited specifier(s)
#   [range|strand|mask_algo_id]
#    * Incompatible with:  entry, range, strand, mask_sequence_with, pig, info,
#   list, recursive, remove_redundant_dbs, list_outfmt,
#   show_blastdb_search_path
# -pig <Integer, >=0>
#   PIG to retrieve
#    * Incompatible with:  entry, entry_batch, target_only, info, list,
#   recursive, remove_redundant_dbs, list_outfmt, show_blastdb_search_path
# -info
#   Print BLAST database information
#    * Incompatible with:  entry, entry_batch, outfmt, strand, target_only,
#   ctrl_a, get_dups, pig, range, mask_sequence, list, remove_redundant_dbs,
#   recursive, list_outfmt, list, recursive, remove_redundant_dbs, list_outfmt,
#   show_blastdb_search_path
#
# *** Sequence retrieval configuration options
# -range <String>
#   Range of sequence to extract in 1-based offsets (Format: start-stop, for
#   start to end of sequence use start - )
#    * Incompatible with:  entry_batch, info, list, recursive,
#   remove_redundant_dbs, list_outfmt, show_blastdb_search_path
# -strand <String, `minus', `plus'>
#   Strand of nucleotide sequence to extract
#   Default = `plus'
#    * Incompatible with:  entry_batch, info, list, recursive,
#   remove_redundant_dbs, list_outfmt, show_blastdb_search_path
# -mask_sequence_with <Integer>
#   Produce lower-case masked FASTA using the algorithm ID specified
#    * Incompatible with:  entry_batch
#
# *** Output configuration options
# -out <File_Out>
#   Output file name
#   Default = `-'
# -outfmt <String>
#   Output format, where the available format specifiers are:
#   		%f means sequence in FASTA format
#   		%s means sequence data (without defline)
#   		%a means accession
#   		%g means gi
#   		%o means ordinal id (OID)
#   		%i means sequence id
#   		%t means sequence title
#   		%l means sequence length
#   		%h means sequence hash value
#   		%T means taxid
#   		%e means membership integer
#   		%L means common taxonomic name
#   		%S means scientific name
#   		%P means PIG
#   		%m means sequence masking data.
#   		   Masking data will be displayed as a series of 'N-M' values
#   		   separated by ';' or the word 'none' if none are available.
#   	If '%f' is specified, all other format specifiers are ignored.
#   	For every format except '%f', each line of output will correspond
#   	to a sequence.
#   Default = `%f'
#    * Incompatible with:  info, list, recursive, remove_redundant_dbs,
#   list_outfmt, show_blastdb_search_path
# -target_only
#   Definition line should contain target entry only
#    * Incompatible with:  pig, info, get_dups, list, recursive,
#   remove_redundant_dbs, list_outfmt, show_blastdb_search_path
# -get_dups
#   Retrieve duplicate accessions
#    * Incompatible with:  info, target_only, list, recursive,
#   remove_redundant_dbs, list_outfmt, show_blastdb_search_path
#
# *** Output configuration options for FASTA format
# -line_length <Integer, >=1>
#   Line length for output
#   Default = `80'
#    * Incompatible with:  list, recursive, remove_redundant_dbs, list_outfmt,
#   show_blastdb_search_path
# -ctrl_a
#   Use Ctrl-A as the non-redundant defline separator
#    * Incompatible with:  info, list, recursive, remove_redundant_dbs,
#   list_outfmt, show_blastdb_search_path
#
# *** BLAST database configuration and discovery options
# -show_blastdb_search_path
#   Displays the default BLAST database search paths
#    * Incompatible with:  entry, entry_batch, outfmt, strand, target_only,
#   ctrl_a, get_dups, pig, range, db, info, mask_sequence, line_length, list,
#   recursive, list_outfmt, remove_redundant_dbs
# -list <String>
#   List BLAST databases in the specified directory
#    * Incompatible with:  info, entry, entry_batch, outfmt, strand,
#   target_only, ctrl_a, get_dups, pig, range, db, info, mask_sequence,
#   line_length, show_blastdb_search_path
# -remove_redundant_dbs
#   Remove the databases that are referenced by another alias file in the
#   directory in question
#    * Incompatible with:  info, entry, entry_batch, outfmt, strand,
#   target_only, ctrl_a, get_dups, pig, range, db, info, mask_sequence,
#   line_length, show_blastdb_search_path
# -recursive
#   Recursively traverse the directory structure to list available BLAST
#   databases
#    * Incompatible with:  info, entry, entry_batch, outfmt, strand,
#   target_only, ctrl_a, get_dups, pig, range, db, info, mask_sequence,
#   line_length, show_blastdb_search_path
# -list_outfmt <String>
#   Output format for the list option, where the available format specifiers
#   are:
#   		%f means the BLAST database absolute file name path
#   		%p means the BLAST database molecule type
#   		%t means the BLAST database title
#   		%d means the date of last update of the BLAST database
#   		%l means the number of bases/residues in the BLAST database
#   		%n means the number of sequences in the BLAST database
#   		%U means the number of bytes used by the BLAST database
#   	For every format each line of output will correspond to a BLAST database.
#   Default = `%f %p'
#    * Incompatible with:  info, entry, entry_batch, outfmt, strand,
#   target_only, ctrl_a, get_dups, pig, range, db, info, mask_sequence,
#   line_length, show_blastdb_search_path
# -logfile <File_Out>
#   File to which the program log should be redirected
