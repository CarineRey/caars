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

from lib import PhyloPrograms
from lib import Aligner

start_time = time.time()

### Option defining
parser = argparse.ArgumentParser(prog = "SeqIntegrator.py",
                                 description='''
    Add sequences to an alignment and merge sequences from defined species if they are phylogenetically enough close.''')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')


##############
requiredOptions = parser.add_argument_group('Required arguments')
requiredOptions.add_argument('-ali', '--alignment', type=str,
                             help='Alignment file name.', required=True)
requiredOptions.add_argument('-fst', '--fasta', type=str,
                             help='Fasta file name', required=True)
requiredOptions.add_argument('-s2t', '--sequence2taxon', type=str,
                             help='Link file name. A tabular file, each line correspond to a sequence name and its species. ', required=True)
requiredOptions.add_argument('-out', '--output_prefix',  type=str, default = "./output",
                   help = "Output prefix (Default ./output)")
##############


##############
Options = parser.add_argument_group('Options')
Options.add_argument('--realign_ali',  action='store_true',
                    help = "A fasta file will be created at each iteration. (default: False)")
Options.add_argument('-tmp',  type=str,
                    help = "Directory to stock all intermediary files for the apytram run. (default: a directory in /tmp which will be removed at the end)",
                    default = "" )
Options.add_argument('-log', type=str, default="seqintegrator.log",
                   help = "a log file to report avancement (default: seq_integrator.log)")

### Option parsing
args = parser.parse_args()

### Read arguments
StartingAlignment = args.alignment
StartingFasta = args.fasta
Seq2Taxon = args.sequence2taxon


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
    TmpDirName = tempfile.mkdtemp(prefix='tmp_SeqIntegrator')

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
if not os.path.isfile(args.alignment):
    logger.error(args.alignment+" (-ali) is not a file.")
    sys.exit(1)
if not os.path.isfile(args.fasta):
    logger.error(args.alignment+" (-fst) is not a file.")
    sys.exit(1)
if not os.path.isfile(args.sequence2taxon):
    logger.error(args.sequence2taxon+" (-s2t) is not a file.")
    sys.exit(1) 
    
    
### Add the fasta file to the existing alignment
logger.info("Add the fasta file to the existing alignment")
MafftProcess = Aligner.Mafft(StartingAlignment)
MafftProcess.AddOption = StartingFasta
MafftProcess.AdjustdirectionOption = False
MafftProcess.AutoOption = True
MafftProcess.OutputFile = "%s/StartMafft.fa" %TmpDirName
if os.path.isfile(StartingAlignment) and os.path.isfile(StartingFasta):
	MafftProcess.get_output()
else:
	logger.error("%s or %s is not a file" %(StartingAlignment,StartingFasta))
	sys.exit(1)

#Remove _R_ add by mafft adjustdirection option
#os.system("sed -i s/_R_//g %s" %MafftProcess.OutputFile)

### Built a tree with the global alignment
logger.info("Built a tree with the global alignment")
FasttreeProcess = PhyloPrograms.Fasttree(MafftProcess.OutputFile)
FasttreeProcess.Nt = True
FasttreeProcess.Gtr = True
FasttreeProcess.OutputTree = "%s/StartTree.tree" %TmpDirName
if os.path.isfile(MafftProcess.OutputFile):
	FasttreeProcess.get_output()
else:
	logger.error("%s is not a file. There was an issue with the previous step." %(MafftProcess.OutputFile))
	sys.exit(1)

### Use phylomerge to merge sequence from a same species
logger.info("Use phylomerge to merge sequence from a same species")
PhylomergeProcess = PhyloPrograms.Phylomerge(MafftProcess.OutputFile, FasttreeProcess.OutputTree)
PhylomergeProcess.SequenceToTaxon = Seq2Taxon
PhylomergeProcess.RearrangeTree = True
PhylomergeProcess.BootstrapThreshold = 0.8
PhylomergeProcess.OutputSequenceFile = "%s/Merged.fa" %TmpDirName
PhylomergeProcess.launch()
if os.path.isfile(MafftProcess.OutputFile) and \
   os.path.isfile(FasttreeProcess.OutputTree) and \
   os.path.isfile(PhylomergeProcess.SequenceToTaxon) :
	PhylomergeProcess.launch()
else:
	logger.error("%s or %s or %s is not a file. There was an issue with the previous step." \
	 %(MafftProcess.OutputFile, FasttreeProcess.OutputTree,PhylomergeProcess.SequenceToTaxon))
	sys.exit(1)

### Realign the merge alignment
logger.info("Realign the merge alignment")
FinalMafftProcess = Aligner.Mafft(PhylomergeProcess.OutputSequenceFile)
FinalMafftProcess.AutoOption = True
FinalMafftProcess.OutputFile = "%s.fa" %OutPrefixName
if os.path.isfile(FinalMafftProcess.InputFile):
	FinalMafftProcess.get_output()
else:
	logger.error("%s is not a file. There was an issue with the previous step." %(FinalMafftProcess.InputFile))
	sys.exit(1)

### Built a tree with the final alignment
logger.info("Built a tree with the final alignment")
FinalFasttreeProcess = PhyloPrograms.Fasttree(MafftProcess.OutputFile)
FinalFasttreeProcess.Nt = True
FinalFasttreeProcess.Gtr = True
FinalFasttreeProcess.OutputTree = "%s.tree" %OutPrefixName

if os.path.isfile(FinalMafftProcess.OutputFile):
	FinalFasttreeProcess.get_output()
else:
	logger.error("%s is not a file. There was an issue with the previous step." %(FinalMafftProcess.OutputFile))
	sys.exit(1)

        
### Remove tempdir if the option --tmp have not been use
if not args.tmp:
    logger.debug("Remove the temporary directory")
    #Remove the temporary directory :
    if "tmp_SeqIntegrator" in TmpDirName:
        shutil.rmtree(TmpDirName)
        
logger.debug("--- %s seconds ---" % (time.time() - start_time))
