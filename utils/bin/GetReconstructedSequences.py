#!/usr/bin/python
# coding: utf-8

# File: GetReconstructedSequences.py
# Created by: Carine Rey
# Created on: September 2016
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

import sys
import os
import glob
import logging
import re


### Set up the logger
# create logger with 'spam_application'
logger = logging.getLogger('GetReconstructedSequences')
logger.setLevel(logging.DEBUG)
# create file handler which logs even debug messages
# create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG) #WARN
# create formatter and add it to the handlers
formatter = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)
# add the handlers to the logger
logger.addHandler(ch)

logger.debug(" ".join(sys.argv))

if len(sys.argv) != 5:
    logger.error("4 arguments are required")
    sys.exit(1)

ali_dir = sys.argv[1]
sp2seq_dir = sys.argv[2]
RefinedSpecies = set(sys.argv[3].split(","))
out_dir = sys.argv[4]


logger.debug(ali_dir)
logger.debug(sp2seq_dir)
logger.debug(" ".join(RefinedSpecies))
logger.debug(out_dir)

### Check input data
if not os.path.isdir(ali_dir):
    logger.error("The alignment directory %s does not exist", ali_dir)
    sys.exit(1)

if not os.path.isdir(sp2seq_dir):
    logger.debug("The sequences2species diectory %s does not exist", sp2seq_dir)
    sys.exit(1)


### Set up the  output directory
if os.path.isdir(out_dir):
    logger.info("The output directory %s exists", out_dir)
else: # we create the directory
    logger.warning("The output directory %s does not exist, it will be created.", out_dir)
    os.makedirs(out_dir)

### Read all files in sp2seq_dir
Seq2Sp_dict = {}

def read_seq2species_file(Seq2Sp_dict, RefinedSpecies, File):
    fam = os.path.basename(File).split('.')[0]
    if os.path.isfile(File):
        f = open(File, "r")
        for line in f:
            line = line.replace("\n", "")
            (sp, seq) = line.split(":")
            if sp in RefinedSpecies:
                Seq2Sp_dict[seq] = (sp, fam)
        f.close()
    return Seq2Sp_dict

for f in glob.glob("%s/*sp2seq.txt" %sp2seq_dir):
    Seq2Sp_dict = read_seq2species_file(Seq2Sp_dict, RefinedSpecies, f)

### Read all alignment in alignment_dir
def read_ali_file(FastaFile, Seq2Sp_dict, AliDict):
    name = ""
    if os.path.isfile(FastaFile):
        f = open(FastaFile, "r")
        for line in f:
            if re.match('[\s\n]', line):
                pass
            elif re.match('>', line):
                name = line[1:-1]
                if not Seq2Sp_dict.has_key(name):
                    name = ""
                else:
                    AliDict.setdefault(name, [])
            elif name:
                AliDict[name].append(line)
            else:
                pass
        f.close()
    return AliDict


def write_seq(AliDict):
    Fasta_File = "%s/amalgam_sequences.fa" %(out_dir)
    string = []
    for (name, seq) in AliDict.items():
        seq = "".join(seq).replace("-", "").replace("\n", "")
        string.extend([">", name, "\n",
                       '\n'.join(seq[i:i+60] for i in range(0, len(seq), 60)),"\n"])
    
    f = open(Fasta_File, "w")
    f.write("".join(string) + "\n")
    f.close()
        
def write_validated_sp2seq(Seq2Sp_dict):
    SeqSpLink_File = "%s/amalgam_sequences.seq2sp2fam.txt" %(out_dir)
    String = []
    sep = "\t"
    for seq in Seq2Sp_dict.keys():
        (sp, fam) = Seq2Sp_dict[seq]
        String.append("%s%s%s%s%s\n" %(seq, sep, sp, sep, fam))

    f = open(SeqSpLink_File, "w")
    f.write("".join(String) + "\n")
    f.close()

AliDict = {}
for f in glob.glob("%s/*" %ali_dir):
    AliDict = read_ali_file(f, Seq2Sp_dict, AliDict)


write_validated_sp2seq(Seq2Sp_dict)
write_seq(AliDict)


sys.exit(0)
