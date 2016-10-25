#!/usr/bin/python
# coding: utf-8

# File: SeqDispatcher.py
# Created by: Carine Rey
# Created on: October 2016
#
#
# Copyright 2016 Carine Rey
# This software is a computer program whose purpose is to assembly
# sequences from RNA-Seq data (paired-end or single-end) using one or
# more reference homologous sequences.
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use,
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info".
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability.
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and,  more generally, to use and operate it in the
# same conditions as regards security.
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.

import os
import sys
import time
import string
import tempfile
import shutil
import logging
import argparse
import subprocess
import re

import pandas

import BlastPlus



start_time = time.time()

### Option defining
parser = argparse.ArgumentParser(prog="CheckFamily.py",
                                 description='''
    Check all sequences of a fasta file are associated with a unique family, else remove them.''')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')


##############
requiredOptions = parser.add_argument_group('Required arguments')
requiredOptions.add_argument('-i', '--input', type=str,
                             help='fasta file name.', required=True)
requiredOptions.add_argument('-t', '--ref_transcriptome', type=str,
                             help='Target fasta file name', required=True)
requiredOptions.add_argument('-f', '--family', type=str,
                             help='family name', required=True)
requiredOptions.add_argument('-d', '--database', type=str,
                             help='''Database prefix name of the ref transcriptome fasta file.
                              If a database with the same name already exists,
                              the existing database will be kept and the database will NOT be rebuilt.
                              (default=: The database will be build in the temporary directory and will be remove at the end.)''',
                              required=False)
requiredOptions.add_argument('-t2f', '--ref_transcriptome2family', type=str,
                             help='Link file name. A tabular file, each line correspond to a sequence name and its family. ', required=True)
requiredOptions.add_argument('-o', '--output', type=str, default="./output.fa",
                   help="Output name (default= ./output.fa)")
##############


##############
Options = parser.add_argument_group('Options')
Options.add_argument('-e', '--evalue', type=float,
                     help="Evalue threshold of the blastn of the queries on the database of the ref transcriptome. (default= 1e-3)",
                     default=1e-3)
Options.add_argument('-tmp', type=str,
                     help="Directory to stock all intermediary files for the job. (default=: a directory in /tmp which will be removed at the end)",
                     default="")
Options.add_argument('-log', type=str, default="checkfamily.log",
                     help="a log file to report avancement (default=: checkfamily.log)")
##############


##############
MiscellaneousOptions = parser.add_argument_group('Miscellaneous options')
MiscellaneousOptions.add_argument('--debug', action='store_true', default=False,
                   help="debug mode, default False")
##############

### Option parsing
args = parser.parse_args()

### Read arguments
FastaFile = args.input
OutputFasta = args.output
TargetFile = args.ref_transcriptome
ExpectedFamily = args.family
Target2FamilyFilename = args.ref_transcriptome2family

Evalue = args.evalue

### Set up the log directory
if args.log:
    LogDirName = os.path.dirname(args.log)
    if not os.path.isdir(LogDirName) and LogDirName:
        os.makedirs(LogDirName)

### Set up the logger
LogFile = args.log
# create logger
logger = logging.getLogger("main")
logger.setLevel(logging.INFO)
# create file handler which logs even debug messages
fh = logging.FileHandler(LogFile)
# create console handler with a higher log level
ch = logging.StreamHandler()
if args.debug:
    logger.setLevel(logging.DEBUG)
    ch.setLevel(logging.DEBUG)
    fh.setLevel(logging.DEBUG)
else:
    ch.setLevel(logging.WARNING)
    fh.setLevel(logging.INFO)

# create formatter and add it to the handlers
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
fh.setFormatter(formatter)
ch.setFormatter(formatter)
# add the handlers to the logger
logger.addHandler(fh)
logger.addHandler(ch)

logger.info(" ".join(sys.argv))

### Set up the working directory
if args.tmp:
    if os.path.isdir(args.tmp):
        logger.info("The temporary directory %s exists", args.tmp)
    else:
        logger.info("The temporary directory %s does not exist, it will be created", args.tmp)
        os.makedirs(args.tmp)
    TmpDirName = args.tmp
else:
    TmpDirName = tempfile.mkdtemp(prefix='checkfamily')

def end(exit_code):
    ### Remove tempdir if the option --tmp have not been use
    if not args.tmp:
        logger.debug("Remove the temporary directory")
        #Remove the temporary directory :
        if "tmp_SeqDispatcher" in TmpDirName:
            shutil.rmtree(TmpDirName)
    sys.exit(exit_code)

### Set up the output directory
if args.output and os.path.dirname(args.output):
    OutDirName = os.path.dirname(args.output)
    OutName = args.output
    if os.path.isdir(OutDirName):
        logger.info("The output directory %s exists", os.path.dirname(args.output))
    elif OutDirName: # if OutDirName is not a empty string we create the directory
        logger.info("The output directory %s does not exist, it will be created", os.path.dirname(args.output))
        os.makedirs(os.path.dirname(args.output))
else:
    logger.error("The output prefix must be defined")
    end(1)

### Check that input files exist
if not os.path.isfile(args.input):
    logger.error(args.input + " (-q) is not a file.")
    end(1)

if not os.path.isfile(args.ref_transcriptome):
    logger.error(args.ref_transcriptome + " (-t) is not a file.")
    end(1)

if not os.path.isfile(args.ref_transcriptome2family):
    logger.error(args.ref_transcriptome2family + " (-t2f) is not a file.")
    end(1)

### Parse input fasta files
## Get query names
logger.info("Get query names")
BashProcess = subprocess.Popen(["grep", "-e", "^>", args.input],
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
OutBashProcess = BashProcess.communicate()
if not  OutBashProcess[1]:
    QueryNames = OutBashProcess[0].strip().replace(">", "").split("\n")
else:
    logger.error(OutBashProcess[1])
    end(1)

## Get ref_transcriptome sequence names
logger.info("Get ref_transcriptome names")
BashProcess = subprocess.Popen(["grep", "-e", "^>",
                                 args.ref_transcriptome],
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
OutBashProcess = BashProcess.communicate()
if not  OutBashProcess[1]:
    TargetNames = OutBashProcess[0].strip().replace(">", "").split("\n")
else:
    logger.error(OutBashProcess[1])
    end(1)

### Parse the ref_transcriptome2family, create dictionnaries target2family and family2target
Target2FamilyTable = pandas.read_csv(Target2FamilyFilename,
                                      sep=None, engine='python',
                                      header=None,
                                      names=["Target", "Family"])

# Check if there are no missing data
MissingTargets = [target for target in TargetNames if not target in Target2FamilyTable.Target.values]
if MissingTargets:
    logger.warning("These targets are not present in the target2family link file, they will be added:\n\t- %s", "\n\t- ".join(MissingTargets))
    MissingData = pandas.DataFrame({"Target" : MissingTargets,
                                     "Family": MissingTargets})
    Target2FamilyTable = pandas.concat([Target2FamilyTable, MissingData], ignore_index=True)


if len(Target2FamilyTable['Target']) != len(Target2FamilyTable['Target'].unique()):
    logger.error("There are not unique target names")
    end(1)

### Check that there is a target database, otherwise build it
logger.info("Check that there is a target database, otherwise build it")
if not args.database:
    DatabaseName = "%s/Target_DB" %TmpDirName
else:
    DatabaseName = args.database
CheckDatabase_BlastdbcmdProcess = BlastPlus.Blastdbcmd(DatabaseName, "", "")
if not CheckDatabase_BlastdbcmdProcess.is_database():
    logger.info("Database %s does not exist", DatabaseName)
    #Build blast formated database from a fasta file
    if not os.path.isfile(TargetFile):
        logger.error("The fasta file (-t) does not exist.")
        end(1)
    if os.path.isdir(os.path.dirname(DatabaseName)):
        logger.info("Database directory exists")
    else:
        logger.info("Database directory does not exist, we create it")
        os.makedirs(os.path.dirname(DatabaseName))
    # database building
    logger.info(DatabaseName + " database building")
    MakeblastdbProcess = BlastPlus.Makeblastdb(TargetFile, DatabaseName)
    (out, err) = MakeblastdbProcess.launch()
    if err:
        end(1)


CheckDatabase_BlastdbcmdProcess = BlastPlus.Blastdbcmd(DatabaseName, "", "")
if not CheckDatabase_BlastdbcmdProcess.is_database():
    logger.error("Problem in the database building")
    logger.info("Database %s does not exist", DatabaseName)
    end(1)
else:
    logger.info("Database %s exists", DatabaseName)

### Blast the query fasta on the target database
logger.info("Blast the query fasta on the target database")
start_blast_time = time.time()
BlastOutputFile = "%s/Queries_Targets.blast" % (TmpDirName)
BlastnProcess = BlastPlus.Blast("blastn", DatabaseName, FastaFile)
BlastnProcess.Evalue = Evalue
BlastnProcess.Task = "blastn"
BlastnProcess.max_target_seqs = 100
BlastnProcess.max_hsps_per_subject = 1
BlastnProcess.OutFormat = "6"
BlastnProcess.Strand="plus"

# Write an empty output file to be sure
OutputFile = open(OutputFasta, "w")
OutputFile.write("")
OutputFile.close()

# Write blast ouptut in BlastOutputFile if the file does not exist
(out, err) = BlastnProcess.launch(BlastOutputFile)
if err:
    end(1)

if not os.stat(BlastOutputFile).st_size:
    logger.info("Blast found no hit")
    end(0)

logger.debug("blast --- %s seconds ---", str(time.time() - start_blast_time))

### Parse blast results
# Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
FieldNames = ["qid", "tid", "id", "alilen", "mis", "gap", "qstart", "qend", "tstart", "tend", "evalue", "score"]
BlastTable = pandas.read_csv(BlastOutputFile, sep=None, engine='python', header=None, names=FieldNames)

# Get Family for each Target:
BlastTableWithFamilies = pandas.merge(BlastTable, Target2FamilyTable, how='left', left_on=['tid'], right_on=['Target'])

### First: Find the best hit for each Query sequences and check family
logger.info("First Step")
RetainedQuery = []

for Query in QueryNames:
    Query = Query.split()[0]
    TmpTable =  BlastTableWithFamilies[BlastTable.qid == Query]
    if not len(TmpTable.index):
        NoHitList.append(Query)
    else:
        TmpBestScore = max(TmpTable.score)

        BestTargetTable = TmpTable[TmpTable.score == TmpBestScore]
        BestTargetTable.is_copy = False

        TmpFamily = BestTargetTable["Family"].unique()
        Family = TmpFamily[0]

        if len(TmpFamily) > 1:
            logger.info("More than one family can be attributed to %s:\n\t- %s\nIt will be discarded.", Query, "\n\t- ".join(TmpFamily))
        elif Family != ExpectedFamily:
            logger.info("Observed family (%s) is different of the expected family (%s). %s will be discarded.", Family, ExpectedFamily, Query)
        else:
            RetainedQuery.extend(BestTargetTable.qid.values)


if not RetainedQuery:
    end(0)

### Second : Build a query database
# Build a query database
logger.info("Build a query database")
QueryDatabaseName = "%s/Query_DB" %TmpDirName
if not os.path.isfile(FastaFile):
    logger.error("The query fasta file (-i) does not exist.")
    end(1)
# database building
logger.info(QueryDatabaseName + " database building")
MakeblastdbProcess = BlastPlus.Makeblastdb(FastaFile, QueryDatabaseName)
ExitCode = MakeblastdbProcess.launch()

CheckDatabase_BlastdbcmdProcess = BlastPlus.Blastdbcmd(QueryDatabaseName, "", "")
if not CheckDatabase_BlastdbcmdProcess.is_database():
    logger.error("Problem in the database building")
    logger.info("Database %s does not exist", QueryDatabaseName)
    end(1)
else:
    logger.info("Database %s exists", QueryDatabaseName)

## Third step: write a fasta which contained all retained sequences
logger.info("Write output files")

#Get retained sequences names
start_blastdbcmd_time = time.time()
TmpFilename = "%s/%s_sequences_names.txt" %(TmpDirName, Family)
TmpFile = open(TmpFilename, "w")
TmpFile.write("\n".join(RetainedQuery))
TmpFile.close()

#Get retained sequences
BlastdbcmdProcess = BlastPlus.Blastdbcmd(QueryDatabaseName, TmpFilename, "")
BlastdbcmdProcess.OutputFile = OutputFasta

(out, err) = BlastdbcmdProcess.launch()

if err:
    end(1)
logger.debug("blastdbcmd --- %s seconds ---", time.time() - start_blastdbcmd_time)

logger.info("--- %s seconds ---", str(time.time() - start_time))

end(0)
