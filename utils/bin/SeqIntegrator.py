#!/usr/bin/python
# coding: utf-8

# File: SeqIntegrator.py
# Created by: Carine Rey
# Created on: March 2016
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
import tempfile
import shutil
import logging
import argparse
import subprocess

import PhyloPrograms
import Aligner

from ete2 import Tree

start_time = time.time()

### Option defining
parser = argparse.ArgumentParser(prog="SeqIntegrator.py",
                                 description='''
    Add sequences to an alignment and merge sequences from defined species if they are phylogenetically enough close.''')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')


##############
requiredOptions = parser.add_argument_group('Required arguments')
requiredOptions.add_argument('-ali', '--alignment', type=str,
                             help='Alignment file name.', required=True)
requiredOptions.add_argument('-fa', '--fasta', type=str,
                             help='Fasta file names delimited by a coma')
requiredOptions.add_argument('-sp2seq', type=str,
                             help='Link file name. A tabular file, each line correspond to a sequence name and its species. File names delimited by comas.', required=True)
requiredOptions.add_argument('-out', '--output_prefix', type=str, default="./output",
                   help="Output prefix (Default ./output)")
##############


##############
Options = parser.add_argument_group('Options')
Options.add_argument('-sptorefine', type=str, default="",
                    help="A list of species names delimited by commas. These species will be concerned by merging. (default: All species will be concerned)")
Options.add_argument('--no_merge', action='store_true', default=False,
                    help="not use phylo merge to merge sequences. (default: False)")
Options.add_argument('--realign_ali', action='store_true', default=False,
                    help="Realign the ali even if no sequences to add. (default: False)")
Options.add_argument('--resolve_polytomy', action='store_true', default=False,
                    help="resolve polytomy. (default: False)")
Options.add_argument('-tmp', type=str,
                    help="Directory to stock all intermediary files for the job. (default: a directory in /tmp which will be removed at the end)",
                    default="")
Options.add_argument('-log', type=str, default="SeqIntegrator.log",
                   help="a log file to report avancement (default: seq_integrator.log)")
Options.add_argument('--debug', action='store_true', default=False,
                   help="debug mode, default False")
### Option parsing
args = parser.parse_args()



### Read arguments
StartingAlignment = args.alignment
SpToRefine = []
if args.sptorefine:
    SpToRefine = set(args.sptorefine.split(","))

if args.fasta:
    FastaFiles = args.fasta.split(",")
else:
    FastaFiles = []
Sp2SeqFiles = args.sp2seq.split(",")

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
fh.setLevel(logging.INFO)
# create console handler with a higher log level
ch = logging.StreamHandler()
if args.debug:
    ch.setLevel(logging.DEBUG)
    fh.setLevel(logging.DEBUG)
    logger.setLevel(logging.DEBUG)
else:
    ch.setLevel(logging.WARNING)
# create formatter and add it to the handlers
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
fh.setFormatter(formatter)
ch.setFormatter(formatter)
# add the handlers to the logger
logger.addHandler(fh)
logger.addHandler(ch)


logger.debug(sys.argv)
### Set up the working directory
if args.tmp:
    if os.path.isdir(args.tmp):
        logger.info("The temporary directory %s exists", args.tmp)
    else:
        logger.info("The temporary directory %s does not exist, it will be created", args.tmp)
        os.makedirs(args.tmp)
    TmpDirName = args.tmp
else:
    TmpDirName = tempfile.mkdtemp(prefix='tmp_SeqIntegrator')

def end(ReturnCode):
    ### Remove tempdir if the option --tmp have not been use
    if not args.tmp:
        logger.debug("Remove the temporary directory")
        #Remove the temporary directory :
        if "tmp_SeqIntegrator" in TmpDirName:
            shutil.rmtree(TmpDirName)
    sys.exit(ReturnCode)

### Set up the output directory
if args.output_prefix:
    OutDirName = os.path.dirname(args.output_prefix)
    OutPrefixName = args.output_prefix
    if os.path.isdir(OutDirName):
        logger.info("The output directory %s exists", os.path.dirname(args.output_prefix))
    elif OutDirName: # if OutDirName is not a empty string we create the directory
        logger.info("The output directory %s does not exist, it will be created", os.path.dirname(args.output_prefix))
        os.makedirs(os.path.dirname(args.output_prefix))
else:
    logger.error("The output prefix must be defined")
    end(1)

### Check that input files exist
if not os.path.isfile(StartingAlignment):
    logger.error(StartingAlignment+" is not a file.")
    end(1)

StartingFastaFiles = []
for f in FastaFiles:
    if os.path.isfile(f) and os.path.getsize(f) > 0:
        StartingFastaFiles.append(f)

StartingSp2SeqFiles = []
for f in Sp2SeqFiles:
    if os.path.isfile(f) and os.path.getsize(f) > 0:
        logger.debug(f)
        StartingSp2SeqFiles.append(f)

### A function to cat input files
def cat(Files, OutputFile):
    (out, err, Output) = ("", "", "")
    command = ["cat"]
    logger.debug(Files)

    if type(Files) == type([]) and len(Files) == 1:
        Output = Files[0]
    else:
        command.extend(Files)
        logger.debug(" ".join(command))
        p = subprocess.Popen(command,
                           stdout=open(OutputFile, 'w'),
                           stderr=subprocess.PIPE)
        (out, err) = p.communicate()
        if err:
            logger.error(err)
        Output = OutputFile
    return (out, err, Output)

def mv(In, Out):
    (out, err) = ("", "")
    command = ["mv", In, Out]

    logger.debug(" ".join(command))
    p = subprocess.Popen(command,
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE)
    (out, err) = p.communicate()
    if err:
        logger.error(err)

    return (out, err)

def cp(In, Out):
    (out, err) = ("", "")
    command = ["cp", In, Out]

    logger.debug(" ".join(command))
    p = subprocess.Popen(command,
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE)
    (out, err) = p.communicate()
    if err:
        logger.error(err)

    return (out, err)

if args.realign_ali:

    ### Realign the input alignment
    InitialMafftProcess = Aligner.Mafft(StartingAlignment)
    InitialMafftProcess.AutoOption = True
    InitialMafftProcess.QuietOption = True
    InitialMafftProcess.OutputFile = "%s/%s.fa" %(TmpDirName, "RealignAli")

    if os.path.isfile(StartingAlignment):
        logger.info("Realign the input alignment")
        _ = InitialMafftProcess.launch()
        StartingAlignment = InitialMafftProcess.OutputFile
    else:
        logger.error("%s is not a file.", StartingAlignment)
        end(1)

### Concate all  sp2seq files
logger.info("Concate all Sp2Seq files")
Sp2Seq = "%s/StartingSp2Seq.txt" %(TmpDirName)
(out, err, Sp2Seq) = cat(StartingSp2SeqFiles, Sp2Seq)

# Check if their are seqeunces to add
if StartingFastaFiles and Sp2SeqFiles:
    logger.info("Sequences to add")
    logger.debug(StartingFastaFiles)
    logger.debug(Sp2SeqFiles)

    ### Concate all fasta files
    StartingFasta = "%s/StartingFasta.fa" %(TmpDirName)
    logger.info("Concate all fasta files")
    (out, err, StartingFasta) = cat(StartingFastaFiles, StartingFasta)


    ### Add the fasta file to the existing alignment
    logger.info("Add the fasta file to the existing alignment")
    MafftProcess = Aligner.Mafft(StartingAlignment)
    MafftProcess.AddOption = StartingFasta
    MafftProcess.AdjustdirectionOption = False
    MafftProcess.AutoOption = True
    MafftProcess.QuietOption = True
    MafftProcess.OutputFile = "%s/StartMafft.fa" %TmpDirName
    if os.path.isfile(StartingAlignment) and os.path.isfile(StartingFasta):
        (out, err) = MafftProcess.launch()
    else:
        logger.error("%s or %s is not a file", StartingAlignment, StartingFasta)
        end(1)

    if args.no_merge:
        logger.info("no_merge=True, sequences will not be merged.")
        LastAli = "%s.fa" %OutPrefixName
        FinalSp2Seq = "%s.sp2seq.txt" %OutPrefixName
        (out, err) = mv(Sp2Seq, FinalSp2Seq)
        (out, err) = mv(MafftProcess.OutputFile, LastAli)
    else:
        ### Built a tree with the global alignment
        logger.info("Built a tree with the global alignment")
        FasttreeProcess = PhyloPrograms.Fasttree(MafftProcess.OutputFile)
        FasttreeProcess.Nt = True
        FasttreeProcess.Gtr = True
        FasttreeProcess.Gamma = True
        FasttreeProcess.OutputTree = "%s/StartTree.tree" %TmpDirName
        if os.path.isfile(MafftProcess.OutputFile):
            FasttreeProcess.get_output()
        else:
            logger.error("%s is not a file. There was an issue with the previous step.", MafftProcess.OutputFile)
            end(1)

        ### Resolve Polytomy
        StartTreeFilename = FasttreeProcess.OutputTree
        if not os.path.isfile(StartTreeFilename):
            logger.error("%s is not a file. There was an issue with the previous step.", StartTreeFilename)
            end(1)
        if args.resolve_polytomy:
            logger.info("Resolve polytomy")
            t = Tree(StartTreeFilename)
            t.resolve_polytomy(recursive=True)
            t.write(format=0, outfile=StartTreeFilename)
        if not os.path.isfile(StartTreeFilename):
            logger.error("%s is not a file. There was an issue with the previous step.", StartTreeFilename)
            end(1)

        ### Use phylomerge to merge sequence from a same species
        logger.info("Use phylomerge to merge sequence from a same species")
        Int1Sp2Seq = "%s/Int1.sp2seq.txt" %TmpDirName
        PhylomergeProcess = PhyloPrograms.Phylomerge(MafftProcess.OutputFile, StartTreeFilename)
        PhylomergeProcess.TaxonToSequence = Sp2Seq
        PhylomergeProcess.RearrangeTree = True
        PhylomergeProcess.BootstrapThreshold = 0.8
        PhylomergeProcess.OutputSequenceFile = "%s/Merged.fa" %TmpDirName
        PhylomergeProcess.OutputTaxonToSequence = Int1Sp2Seq
        if SpToRefine:
            logger.debug("Species to refine:\n"+"\n".join(SpToRefine))
            SpToRefineFilename = "%s/SpToRefine.txt" %TmpDirName
            SpToRefineFile = open(SpToRefineFilename, "w")
            SpToRefineFile.write("\n".join(SpToRefine)+"\n")
            SpToRefineFile.close()
            PhylomergeProcess.TaxonsToRefine = SpToRefineFilename

        if os.path.isfile(MafftProcess.OutputFile) and \
           os.path.isfile(StartTreeFilename) and \
           os.path.isfile(PhylomergeProcess.TaxonToSequence):
            PhylomergeProcess.launch()
        else:
            logger.error("%s or %s or %s is not a file. There was an issue with the previous step.",
            MafftProcess.OutputFile, StartTreeFilename, PhylomergeProcess.TaxonToSequence)
            end(1)

        logger.info("Realign the merge alignment")
        ### Realign the merged alignment
        Int1MafftProcess = Aligner.Mafft(PhylomergeProcess.OutputSequenceFile)
        Int1MafftProcess.AutoOption = True
        Int1MafftProcess.QuietOption = True
        Int1MafftProcess.OutputFile = "%s/Int1.fa" %TmpDirName
        if os.path.isfile(Int1MafftProcess.InputFile):
            (out, err) = Int1MafftProcess.launch()
        else:
            logger.error("%s is not a file. There was an issue with the previous step. (Phylomerge)", Int1MafftProcess.InputFile)
            end(1)
        Int1Ali = Int1MafftProcess.OutputFile

        ### Built a tree with the realigned merged alignment"
        logger.info("Built a tree with the realigned merged alignment")
        Int1FasttreeProcess = PhyloPrograms.Fasttree(Int1Ali)
        Int1FasttreeProcess.Nt = True
        Int1FasttreeProcess.Gtr = True
        Int1FasttreeProcess.Gamma = True
        Int1FasttreeProcess.OutputTree = "%s/Int1.tree" %TmpDirName
        if os.path.isfile(Int1Ali):
            Int1FasttreeProcess.get_output()
        else:
            logger.error("%s is not a file. There was an issue with the previous step.", Int1Ali)
            end(1)

        ### Resolve Polytomy
        Int1TreeFilename = Int1FasttreeProcess.OutputTree
        if not os.path.isfile(Int1TreeFilename):
            logger.error("%s is not a file. There was an issue with the previous step.", Int1TreeFilename)
            end(1)
        if args.resolve_polytomy:
            logger.info("Resolve polytomy")
            t = Tree(Int1TreeFilename)
            t.resolve_polytomy(recursive=True)
            t.write(format=0, outfile=Int1TreeFilename)
            if not os.path.isfile(Int1TreeFilename):
                logger.error("%s is not a file. There was an issue with the previous step. (Resolve polytomy)", Int1TreeFilename)
                end(1)

        ### Use phylomerge to merge sequence from a same species
        logger.info("Use phylomerge to merge sequence from a same species")
        FinalSp2Seq = "%s.sp2seq.txt" %OutPrefixName
        PhylomergeProcess = PhyloPrograms.Phylomerge(Int1Ali, Int1TreeFilename)
        PhylomergeProcess.TaxonToSequence = Int1Sp2Seq
        PhylomergeProcess.RearrangeTree = True
        PhylomergeProcess.BootstrapThreshold = 0.8
        PhylomergeProcess.OutputSequenceFile = "%s/Int1Merged.fa" %TmpDirName
        PhylomergeProcess.OutputTaxonToSequence = FinalSp2Seq
        if SpToRefine:
            logger.debug("Species to refine:\n"+"\n".join(SpToRefine))
            SpToRefineFilename = "%s/SpToRefine.txt" %TmpDirName
            SpToRefineFile = open(SpToRefineFilename, "w")
            SpToRefineFile.write("\n".join(SpToRefine)+"\n")
            SpToRefineFile.close()
            PhylomergeProcess.TaxonsToRefine = SpToRefineFilename

        if os.path.isfile(Int1Ali) and \
           os.path.isfile(Int1TreeFilename) and \
           os.path.isfile(Int1Sp2Seq):
            PhylomergeProcess.launch()
        else:
            logger.error("%s or %s or %s is not a file. There was an issue with the previous step.",
            Int1Ali, Int1TreeFilename, Int1Sp2Seq)
            end(1)

        logger.info("Realign the merge alignment")
        ### Realign the last alignment
        FinalMafftProcess = Aligner.Mafft(PhylomergeProcess.OutputSequenceFile)
        FinalMafftProcess.AutoOption = True
        FinalMafftProcess.QuietOption = True
        FinalMafftProcess.OutputFile = "%s.fa" %OutPrefixName
        if os.path.isfile(FinalMafftProcess.InputFile):
            (out, err) = FinalMafftProcess.launch()
        else:
            logger.error("%s is not a file. There was an issue with the previous step. (Phylomerge 2nd)", FinalMafftProcess.InputFile)
            end(1)
        LastAli = FinalMafftProcess.OutputFile


else: #No sequences to add
    logger.warning("No sequences to add, the input file will be the output file")
    LastAli = "%s.fa" %OutPrefixName
    FinalSp2Seq = "%s.sp2seq.txt" %OutPrefixName
    (out, err) = cp(Sp2Seq, FinalSp2Seq)
    if args.realign_ali:
        (out, err) = mv(StartingAlignment, LastAli)
    else:
        (out, err) = cp(StartingAlignment, LastAli)

### Built a tree with the final alignment
logger.info("Built a tree with the final alignment")
FinalFasttreeProcess = PhyloPrograms.Fasttree(LastAli)
FinalFasttreeProcess.Nt = True
FinalFasttreeProcess.Gtr = True
FinalFasttreeProcess.Gamma = True
FinalFasttreeProcess.OutputTree = "%s.tree" %OutPrefixName

if os.path.isfile(LastAli):
    FinalFasttreeProcess.get_output()
else:
    logger.error("%s is not a file. There was an issue with the previous step.", LastAli)
    end(1)

### Resolve Polytomy
FinalTreeFilename = FinalFasttreeProcess.OutputTree
if not os.path.isfile(FinalTreeFilename):
    logger.error("%s is not a file. There was an issue with the previous step.", FinalTreeFilename)
    end(1)
if not os.path.getsize(FinalTreeFilename):
    logger.error("%s is empty. There was an issue with the previous step.", FinalTreeFilename)
    end(1)

if args.resolve_polytomy:
    logger.info("Resolve polytomyi in %s", FinalTreeFilename)
    t = Tree(FinalTreeFilename)
    t.resolve_polytomy(recursive=True)
    t.write(format=0, outfile=FinalTreeFilename)


if not os.path.isfile(FinalTreeFilename):
    FinalFasttreeProcess.get_output()
    logger.error("%s is not a file. There was an issue with the previous step.", FinalTreeFilename)
    end(1)

if not os.path.getsize(FinalTreeFilename):
    logger.error("%s is empty. There was an issue with the previous step.", FinalTreeFilename)
    end(1)

logger.debug("--- %s seconds ---", str(time.time() - start_time))
end(0)
