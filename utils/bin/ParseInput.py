#!/usr/bin/python
# coding: utf-8

# File: ParseInput.py
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



#python bin/ParseInput.py Bistro/example_2/Config_files/Config_RNA-seq.tsv Bistro/example_2/Config_files/Species_tree.nw  Bistro/example_2/Alignment_data/ Bistro/example_2/Sequences_Species_links/ Bistro/example_2/

import sys
import os
import glob
import logging
import re

import ete2


### Set up the logger
# create logger with 'spam_application'
logger = logging.getLogger('ParseInput')
logger.setLevel(logging.WARNING)
# create file handler which logs even debug messages
# create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.WARNING) #WARN
# create formatter and add it to the handlers
formatter = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)
# add the handlers to the logger
logger.addHandler(ch)

logger.debug(" ".join(sys.argv))

if len(sys.argv) != 6:
    logger.error("5 arguments are required")
    sys.exit(1)

config_file = sys.argv[1]
species_tree_file = sys.argv[2]
ali_dir = sys.argv[3]
seq2sp_dir = sys.argv[4]
out_dir = sys.argv[5]

logger.debug(config_file)
logger.debug(ali_dir)
logger.debug(seq2sp_dir)
logger.debug(out_dir)

### Check input data
if not os.path.isfile(config_file):
    logger.error("The config file %s does not exist", config_file)
    sys.exit(1)

if not os.path.isfile(species_tree_file):
    logger.error("The species tree file %s does not exist", species_tree_file)
    sys.exit(1)

if not os.path.isdir(ali_dir):
    logger.error("The alignment directory %s does not exist", ali_dir)
    sys.exit(1)

if not os.path.isdir(seq2sp_dir):
    logger.debug("The sequences2species diectory %s does not exist", seq2sp_dir)
    sys.exit(1)


### Set up the  output directory
if os.path.isdir(out_dir):
    logger.info("The output directory %s exists", out_dir)
else: # we create the directory
    logger.warning("The output directory %s does not exist, it will be created.", out_dir)
    os.makedirs(out_dir)


### Read the species tree
t = ete2.Tree(species_tree_file)
All_Species = [leaf.name for leaf in t.iter_leaves()]

logger.info("Sp : %s", All_Species)

### Retrieve Reference species for Trinity or apytram:
RefSpTrinity = []
RefSpApytram = []
RnaSp = []
f = open(config_file, "r")
HeaderConf = f.readline()
for line in f:
    line_list = line.split("\t")
    if len(line_list) == 10:
        (rna_id, sp, ref_species, path_fastq_single, path_fastq_left, path_fastq_right, orientation, run_trinity, path_assembly, run_apytram) = line_list
        if path_assembly != "-":
            if not os.path.isfile(path_assembly):
                logger.error("The given trinity assembly file %s does not exist for %s", path_assembly, rna_id)
                sys.exit(1)
        if sp in All_Species and ref_species in All_Species:
            if run_apytram.strip() in ["y", "yes", "Y", "Yes"]:
                RefSpApytram.extend(ref_species.split(","))
            if run_trinity.strip() in ["y", "yes", "Y", "Yes"]:
                RefSpTrinity.extend(ref_species.split(","))
            if sp not in RnaSp:
                RnaSp.append(sp)
            if orientation not in ["FR", "RF", "F", "R", "US", "UP"]:
                logger.error("orientation must be  in [FR,RF,F,R,US,UP] and not: %s", orientation)
                sys.exit(1)
        else:
            logger.error("%s or %s are not in the species tree.\nSpecies in the species tree:%s", sp, ref_species, "\n\t".join(All_Species))
            sys.exit(1)

    else:
        logger.error("Config file has not 10 elements in line:\n%s", line)
        sys.exit(1)
logger.info("Ref Trinity : %s", RefSpTrinity)
logger.info("Ref Apytram : %s", RefSpApytram)
logger.info("Rna species : %s", RnaSp)

### Read all files in seq2sp_dir
Seq2Sp_dict = {}

def read_seq2species_file(Seq2Sp_dict, File):
    if os.path.isfile(File):
        f = open(File, "r")
        for line in f:
            (seq, sp) = line.split("\t")
            if not Seq2Sp_dict.has_key(seq):
                Seq2Sp_dict[seq] = sp.replace("\n", "")
            else:
                logger.debug("ERROR : Sequence name \"%s\" is not unique")
                sys.exit(1)
        f.close()
    return Seq2Sp_dict

for f in glob.glob("%s/*.tsv" %seq2sp_dir):
    Seq2Sp_dict = read_seq2species_file(Seq2Sp_dict, f)

if len(set(Seq2Sp_dict.values())) == len(set(Seq2Sp_dict.values()).intersection(set(All_Species))):
    logger.info("All species in the species tree have at least one sequence in a Seq2Sp file")
else:
    logger.error("There is not the same number of species in the species tree and in the Seq2sp directory")
    sys.exit(1)


### Read all alignment in alignment_dir

def read_ali_file(FastaFile):
    AliDict = {}
    name = 0
    err = 1
    if os.path.isfile(FastaFile):
        err = 0
        f = open(FastaFile, "r")
        for line in f:
            if re.match('[\s\n]', line):
                pass
            elif re.match('>', line):
                name = line[1:-1]
            else:
                AliDict.setdefault(name, []).append(line)
        f.close()
    if AliDict.has_key(0):
        err = 1
    return AliDict, err


# Write a transcriptome for each RefSpTrinity
TranscriptomePath_dict = {}
TranscriptomeDirPath = "%s/R_Sp_transcriptomes" %(out_dir)
if not os.path.isdir(TranscriptomeDirPath):
    os.makedirs(TranscriptomeDirPath)

SeqSpLinkDirPath = "%s/Alignments_Species2Sequences" %(out_dir)
if not os.path.isdir(SeqSpLinkDirPath):
    os.makedirs(SeqSpLinkDirPath)

SeqFamLinkDirPath = "%s/R_Sp_Seq_Fam_links" %(out_dir)
if not os.path.isdir(SeqFamLinkDirPath):
    os.makedirs(SeqFamLinkDirPath)

ApytramGeneFamDirPath = "%s/R_Sp_Gene_Families" %(out_dir)
if not os.path.isdir(ApytramGeneFamDirPath):
    os.makedirs(ApytramGeneFamDirPath)


def write_validated_sp2seq(SeenSeq2SpDict_i, Family):
    SeqSpLink_File = "%s/alignments.%s.sp2seq.txt" %(SeqSpLinkDirPath, Family)
    String = []
    sep = ":"
    for seq in SeenSeq2SpDict_i.keys():
        String.append("%s%s%s\n" %(SeenSeq2SpDict_i[seq], sep, seq))

    f = open(SeqSpLink_File, "w")
    f.write("".join(String))
    f.close()


def write_seq_ref_Trinity(Ref_dic_trinity, AliDict_i, Family):
    for sp in Ref_dic_trinity.keys():
        #Transcriptome:
        Transcriptome_File = "%s/%s_transcriptome.fa" %(TranscriptomeDirPath, sp)
        string = []
        for name in Ref_dic_trinity[sp]:
            seq = ''.join(AliDict_i[name]).replace("-", "")
            string.extend([">", name, "\n",
                           '\n'.join(seq[i:i+60] for i in range(0, len(seq), 60))])

        f = open(Transcriptome_File, "a")
        f.write("".join(string))
        f.close()

        #Tab Seq 2 Fam:
        FamSeqLink_File = "%s/%s_Fam_Seq.fa" %(SeqFamLinkDirPath, sp)
        string = []
        for name in Ref_dic_trinity[sp]:
            string.extend([name, "\t", Family, "\n"])

        f = open(FamSeqLink_File, "a")
        f.write("".join(string))
        f.close()

def write_seq_ref_apytram(Ref_dic_trinity, AliDict_i, Family):
    for sp in Ref_dic_apytram.keys():
        #gene family:
        GeneFamily_File = "%s/%s.%s.fa" %(ApytramGeneFamDirPath, sp, Family)
        string = []
        for name in Ref_dic_trinity[sp]:
            seq = ''.join(AliDict_i[name]).replace("-", "")
            string.extend([">", name, "\n",
                           '\n'.join(seq[i:i+60] for i in range(0, len(seq), 60))])

        f = open(GeneFamily_File, "w")
        f.write("".join(string))
        f.close()

SeenSeq2SpDict = {}
CountDict2 = {}
Nb_Family = 0
for f in glob.glob("%s/*" %ali_dir):
    Family = os.path.basename(f).split('.')[0]
    Nb_Family += 1
    AliDict_i, err = read_ali_file(f)
    if err:
        logger.error("%s is not a fasta file", f)
        sys.exit(1)
    # Check all sequence name in Seq2SpDict
    if not len(AliDict_i.keys()) == len(set(AliDict_i.keys()).intersection(set(Seq2Sp_dict.keys()))):
        logger.error("All sequences present in %f are not in a file from %s", f, seq2sp_dir)
        sys.exit(1)
    # Check all sp in All_species and write each temporary files
    Ref_dic_trinity = {}
    Ref_dic_apytram = dict([(key, []) for key in RefSpApytram])
    SeenSeq2SpDict_i = {}

    for seq in AliDict_i.keys():
        sp = Seq2Sp_dict[seq]
        if SeenSeq2SpDict.has_key(seq):
            logger.error("Sequence name:%s is not unique")
            sys.exit(1)
        SeenSeq2SpDict[seq] = sp
        SeenSeq2SpDict_i[seq] = sp
        CountDict2.setdefault(sp, {"Nb_seq":0, "Nb_family":0, "Families":[]})
        CountDict2[sp]["Nb_seq"] += 1
        if not Family in CountDict2[sp]["Families"]:
            CountDict2[sp]["Nb_family"] += 1
            CountDict2[sp]["Families"].append(Family)
        if sp in RefSpTrinity + RefSpApytram:
            if sp in RefSpTrinity:
                Ref_dic_trinity.setdefault(sp, []).append(seq)
            if sp in RefSpApytram:
                Ref_dic_apytram.setdefault(sp, []).append(seq)
        elif not sp in All_Species:
            logger.error("%s not in the species tree (%s)", Seq2Sp_dict[seq], species_tree_file)
            sys.exit(1)

    write_validated_sp2seq(SeenSeq2SpDict_i, Family)
    write_seq_ref_Trinity(Ref_dic_trinity, AliDict_i, Family)
    write_seq_ref_apytram(Ref_dic_apytram, AliDict_i, Family)


# Statistics:
# Number of sequences by species:
logger.info("Number of sequence by species:\n"+
            "\n".join(["%s: %s" %(s, c["Nb_seq"]) for (s, c) in CountDict2.items()]))

logger.info("Number of family by species:\n"+
            "\n".join(["%s: %s" %(s, c["Nb_family"]) for (s, c) in CountDict2.items()]))
logger.info("Number of families: %s", Nb_Family)

for sp in RefSpTrinity + RefSpApytram:
    if not CountDict2[sp]:
        logger.error("No sequence for  %s", sp)
        sys.exit(1)

sys.exit(0)
