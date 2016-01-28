#!/usr/bin/python
# coding: utf-8
import os
import sys
import time
import tempfile
import shutil
import logging
import argparse
import subprocess
import re

import pandas

from lib import PhyloPrograms
from lib import Aligner
from lib import BlastPlus



start_time = time.time()

### Option defining
parser = argparse.ArgumentParser(prog = "SeqDispatcher.py",
                                 description='''
    Attribute a family to a list of sequences according to a given target
    sequences list wich have an asocciated family. ''')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')


##############
requiredOptions = parser.add_argument_group('Required arguments')
requiredOptions.add_argument('-q', '--query', type=str,
                             help='Query fasta file name.', required=True)
requiredOptions.add_argument('-t', '--target', type=str,
                             help='Target fasta file name', required=True)
requiredOptions.add_argument('-d', '--database', type=str,
                             help='''Database prefix name of the target fasta file.
                              If a database with the same name already exists,
                              the existing database will be kept and the database will NOT be rebuilt.
                              (Default: The database will be buil in the temporary directory and will be remove at the end.)''',
                              required=False)
requiredOptions.add_argument('-t2f', '--target2family', type=str,
                             help='Link file name. A tabular file, each line correspond to a sequence name and its family. ', required=True)
requiredOptions.add_argument('-out', '--output_prefix',  type=str, default = "./output",
                   help = "Output prefix (Default ./output)")
##############


##############
Options = parser.add_argument_group('Options')
Options.add_argument('-e', '--evalue',  type=float,
                    help = "Evalue threshold of the blastn of the queries on the database of targets. (Default 1e-3)",
                    default = 1e-3 )

Options.add_argument('-tmp',  type=str,
                    help = "Directory to stock all intermediary files for the apytram run. (default: a directory in /tmp which will be removed at the end)",
                    default = "" )
Options.add_argument('-log', type=str, default="seq_dispatcher.log",
                   help = "a log file to report avancement (default: seq_dispatcher.log)")
##############


##############
MiscellaneousOptions = parser.add_argument_group('Miscellaneous options')
MiscellaneousOptions.add_argument('-threads',  type=int,
                    help = "Number of available threads. (Default 1)",
                    default = 1 )
##############

### Option parsing
args = parser.parse_args()

### Read arguments
QueryFile = args.query
TargetFile = args.target
Target2FamilyFilename = args.target2family

Evalue = args.evalue
Threads = args.threads

### Set up the log directory
if args.log:
    LogDirName = os.path.dirname(args.log)
    if not os.path.isdir(LogDirName) and LogDirName:
        os.makedirs(LogDirName)

### Set up the logger
LogFile = args.log
# create logger
logger = logging.getLogger("main")
logger.setLevel(logging.DEBUG)
# create file handler which logs even debug messages
fh = logging.FileHandler(LogFile)
fh.setLevel(logging.DEBUG)
# create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
# create formatter and add it to the handlers
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
fh.setFormatter(formatter)
ch.setFormatter(formatter)
# add the handlers to the logger
logger.addHandler(fh)
logger.addHandler(ch)

### Set up the working directory
if args.tmp:
    if os.path.isdir(args.tmp):
        logger.info("The temporary directory %s exists" %(args.tmp) )
    else:
        logger.info("The temporary directory %s does not exist, it will be created" % (args.tmp))
        os.makedirs(args.tmp)
    TmpDirName = args.tmp
else:
    TmpDirName = tempfile.mkdtemp(prefix='tmp_SeqDispatcher')

### Set up the output directory
if args.output_prefix:
    OutDirName = os.path.dirname(args.output_prefix)
    OutPrefixName = args.output_prefix
    if os.path.isdir(OutDirName):
        logger.info("The output directory %s exists" %(os.path.dirname(args.output_prefix)) )
    elif OutDirName: # if OutDirName is not a empty string we create the directory
        logger.info("The temporary directory %s does not exist, it will be created" % (os.path.dirname(args.output_prefix)))
        os.makedirs(os.path.dirname(args.output_prefix))
else:
    logger.error("The output prefix must be defined")
    sys.exit(1)

### Check that input files exist
if not os.path.isfile(args.query):
    logger.error(args.alignment+" (-q) is not a file.")
    sys.exit(1)

if not os.path.isfile(args.target):
    logger.error(args.target+" (-t) is not a file.")
    sys.exit(1)
        
if not os.path.isfile(args.target2family):
    logger.error(args.target2family+" (-t2f) is not a file.")
    sys.exit(1) 

### Parse input fasta files
## Get query names
logger.info("Get query names")
BashProcess = subprocess.Popen(["grep","-e", "^>",args.query],
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
OutBashProcess = BashProcess.communicate()
if not  OutBashProcess[1]:
    QueryNames = OutBashProcess[0].strip().replace(">","").split("\n")
else:
    logger.error(OutBashProcess[1])
    sys.exit(1)

## Get target names
logger.info("Get target names")
BashProcess = subprocess.Popen(["grep","-e", "^>",args.target],
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
OutBashProcess = BashProcess.communicate()
if not  OutBashProcess[1]:
    TargetNames = OutBashProcess[0].strip().replace(">","").split("\n")
else:
    logger.error(OutBashProcess[1])
    sys.exit(1)
                
### Parse the target2family, create dictionnaries target2family and family2target
Target2FamilyDic = {}
Family2TargetDic = {}

Target2FamilyTable = pandas.read_csv(Target2FamilyFilename, sep = None , engine='python', header = None, names = ["Target","Family"])

# Check if there are no missing data
MissingTargets = [target for target in TargetNames if not target in Target2FamilyTable.Target.values]
if MissingTargets:
    logger.warning("These targets are not present in the target2family link file, they will be added:\n\t- %s" %("\n\t- ".join(MissingTargets)))
    MissingData = pandas.DataFrame({"Target" : MissingTargets,
                                     "Family": MissingTargets})
    Target2FamilyTable = pandas.concat([Target2FamilyTable,MissingData],ignore_index=True)

    
if len(Target2FamilyTable['Target']) != len(Target2FamilyTable['Target'].unique()):
    logger.error("There are not unique target names")
    sys.exit(1)
Target2FamilyDic = Target2FamilyTable.set_index('Target').T.to_dict('list')
    
    
for family in Target2FamilyTable['Family'].unique():
    Family2TargetDic[family] = Target2FamilyTable['Target'][Target2FamilyTable['Family'] == family].values
    #print "Family: %s\tNumber of targets: %s" %(family, len(Target2FamilyTable.Target[Target2FamilyTable['Family'] == family]))


### Check that there is a target database, otherwise build it
logger.info("Check that there is a target database, otherwise build it")
if not args.database:
    DatabaseName = "%s/Target_DB" %TmpDirName
else:
    DatabaseName = args.database
CheckDatabase_BlastdbcmdProcess = BlastPlus.Blastdbcmd(DatabaseName, "")
if not CheckDatabase_BlastdbcmdProcess.is_database():
    logger.info("Database %s does not exist" % DatabaseName)
    #Build blast formated database from a fasta file
    if not os.path.isfile(TargetFile):
        logger.error("The fasta file (-t) does not exist.")
        sys.exit(1)
    if os.path.isdir(os.path.dirname(DatabaseName)):
        logger.info("Database directory exists")
    else:
        logger.info("Database directory does not exist, we create it")
        os.makedirs(os.path.dirname(DatabaseName))
    # database building
    logger.info(DatabaseName + " database building")
    MakeblastdbProcess = BlastPlus.Makeblastdb(TargetFile,DatabaseName)
    ExitCode = MakeblastdbProcess.launch()

CheckDatabase_BlastdbcmdProcess = BlastPlus.Blastdbcmd(DatabaseName, "")
if not CheckDatabase_BlastdbcmdProcess.is_database():
    logger.error("Problem in the database building")
    logger.info("Database %s does not exist" % DatabaseName)
    sys.exit(1)
else:
    logger.info("Database %s exists" % DatabaseName)

### Blast the query fasta on the target database
logger.info("Blast the query fasta on the target database")
start_blast_time = time.time()
BlastOutputFile = "%s/Queries_Targets.blast" % (TmpDirName)
BlastnProcess = BlastPlus.Blast("blastn", DatabaseName, QueryFile)
BlastnProcess.Evalue = Evalue
BlastnProcess.Task = "megablast"
BlastnProcess.max_target_seqs = 500
BlastnProcess.Threads = Threads
BlastnProcess.OutFormat = "6"
# Write blast ouptut in BlastOutputFile if the file does not exist
if not os.path.isfile(BlastOutputFile):
    ExitCode = BlastnProcess.launch(BlastOutputFile)
else:
    logger.warn("%s has already been created, it will be used" %BlastOutputFile )
logger.debug("blast --- %s seconds ---" % (time.time() - start_blast_time))

### Parse blast results
# Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
FieldNames = ["qid", "tid", "id", "alilen", "mis", "gap", "qstart", "qend", "tstart", "tend", "evalue", "score"]
BlastTable = pandas.read_csv(BlastOutputFile, sep = None , engine='python', header = None, names = FieldNames)

# First: Find the best hit for each Query sequences and create a Hit dictonary
logger.info("First Step")
HitDic = {}
NoHitList = []
for Query in QueryNames:
    Query = Query.split()[0]
    TmpTable = BlastTable[BlastTable.qid == Query]
    if not len(TmpTable.index):
        NoHitList.append(Query)
    else:
        TmpBestScore = max(TmpTable.score)
        BestTarget = TmpTable.tid[TmpTable.score == TmpBestScore].values
        TmpFamily = []

        for Target in BestTarget:
            Family = Target2FamilyDic[Target][0]
            if not Family in TmpFamily:
                TmpFamily.append(Family)
                
            if not HitDic.has_key(Family):
                HitDic[Family] = {}
                
            if HitDic[Family].has_key(Target):
                HitDic[Family][Target]["Query"].append(Query)
                HitDic[Family][Target]["Score"].append(TmpBestScore)
            else:
                HitDic[Family][Target] = {"Query":[Query],"Score":[TmpBestScore]}
        
        if len(TmpFamily) > 1:
            logger.warning("More than one family has been attributed to %s:\n\t- %s" %(Query,"\n\t- ".join(TmpFamily)))
       
if NoHitList:
    logger.debug("Queries wihout blast hit:\n\t- %s" %("\n\t- ".join(NoHitList)))                

# Second: For each family, for each target with an hit we kept hits with a score >=0.9 of the best hit
logger.info("Second Step")
for Family in HitDic.keys():
    for Target in HitDic[Family]:
        BestScore = max(HitDic[Family][Target]["Score"])
        L = len(HitDic[Family][Target]["Score"])
        Threshold = 0.9
        HitDic[Family][Target]["Retained"] = [HitDic[Family][Target]["Query"][i] for i in range(L) if  HitDic[Family][Target]["Score"][i] >= (Threshold*BestScore)]
            
### Write output file
## Build a query database
logger.info("Build a query database")
QueryDatabaseName = "%s/Query_DB" %TmpDirName
if not os.path.isfile(QueryFile):
    logger.error("The query fasta file (-q) does not exist.")
    sys.exit(1)
# database building
logger.info(QueryDatabaseName + " database building")
MakeblastdbProcess = BlastPlus.Makeblastdb(QueryFile,QueryDatabaseName)
ExitCode = MakeblastdbProcess.launch()

CheckDatabase_BlastdbcmdProcess = BlastPlus.Blastdbcmd(QueryDatabaseName, "")
if not CheckDatabase_BlastdbcmdProcess.is_database():
    logger.error("Problem in the database building")
    logger.info("Database %s does not exist" % QueryDatabaseName)
    sys.exit(1)
else:
    logger.info("Database %s exists" % QueryDatabaseName)


## Third step: For each family, write a fasta which contained all retained family
logger.info("Write output files")
OutputTableString = ""
for Family in HitDic.keys():
    logger.debug("for %s" %Family)
    TmpRetainedNames = []
    TmpDic = {}
    for Target in HitDic[Family]:
        TmpRetainedNames.extend(HitDic[Family][Target]["Retained"])
        for Query in HitDic[Family][Target]["Retained"]:
            TmpDic[Query] = Target
    ## Write outputs
    #Get retained sequences names
    start_blastdbcmd_time = time.time()
    TmpFilename = "%s/%s_sequences_names.txt" %(TmpDirName, Family)
    TmpFile =open(TmpFilename,"w")
    TmpFile.write("\n".join(TmpRetainedNames))
    TmpFile.close()
    #Get retained sequences
    FamilyOutputName = "%s_%s.fasta" %(OutPrefixName, Family)
    BlastdbcmdProcess = BlastPlus.Blastdbcmd(QueryDatabaseName, TmpFilename)
    FastaString = BlastdbcmdProcess.get_output()
    print FastaString
    logger.debug("blastdbcmd --- %s seconds ---" %(time.time() - start_blastdbcmd_time))
    #Rename sequences:
    NewString = ""
    for line in FastaString.strip().split("\n"):
        if re.match(">",line):
            OldName = line
            Name = line.replace(">lcl|","").strip().split()[0]
            Target = TmpDic[Name]
            NewName = ">%s|%s" %(Name,Target)
            NewString += NewName+"\n"
            # Write in the output table query target family
            OutputTableString += "%s\t%s\t%s\n" %(Name,Target,Family)
        else:
            NewString += line +"\n"
    
    FamilyOutput = open(FamilyOutputName,"w")
    FamilyOutput.write(NewString)
    FamilyOutput.close()

OutputTableFilename = "%s_table.tsv" %(OutPrefixName, Family)
OutputTableFile = open(OutputTableFilename,"w")
OutputTableFile.write(OutputTableString)
OutputTableFile.close()
        
### Remove tempdir if the option --tmp have not been use
if not args.tmp:
    logger.debug("Remove the temporary directory")
    #Remove the temporary directory :
    if "tmp_apytram" in TmpDirName:
        shutil.rmtree(TmpDirName)

logger.info("--- %s seconds ---" % (time.time() - start_time))
