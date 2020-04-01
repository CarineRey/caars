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

import ete3

sep = ":"

### Set up the logger
# create logger with 'spam_application'
logger = logging.getLogger('ParseInput')
logger.setLevel(logging.WARN)
# create file handler which logs even debug messages
# create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.WARN) #WARN
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
    logger.error("The sequences2species diectory %s does not exist", seq2sp_dir)
    sys.exit(1)


### Set up the  output directory
if os.path.isdir(out_dir):
    logger.debug("The output directory %s exists", out_dir)
else: # we create the directory
    logger.warning("The output directory %s does not exist, it will be created.", out_dir)
    os.makedirs(out_dir)


### Read the species tree
t = ete3.Tree(species_tree_file)
All_Species = [leaf.name for leaf in t.iter_leaves()]

for sp in All_Species:
    if sp.replace(sep, "") != sp :
        logger.error("The species [%s] contained the separator [%s] used in CAARS. Please remove it." %(sp, sep))
        sys.exit(1)

logger.info("Sp:\n%s", ";".join(All_Species))

### Retrieve Reference species for Trinity or apytram:
RefSpTrinity = []
RefSpApytram = []
RnaSp = []

error_nb = 0
l_nb = 0
logger.info("Parse the sample sheet")
with open(config_file, "r") as f:
    HeaderConf = f.readline()
    for line in f:
        l_nb +=1
        line_list = line.strip().split("\t")
        if line_list == [""]:
            logger.warning("l%i: empty line", l_nb)
        elif len(line_list) == 11:
            (rna_id, sp, apytram_group, ref_species, path_fastq_single, path_fastq_left, path_fastq_right, orientation, run_trinity, path_assembly, run_apytram) = line_list
            if path_assembly != "-":
                if not os.path.isfile(path_assembly):
                    error_nb += 1
                    logger.error("l%i: The given trinity assembly file %s does not exist for %s", l_nb, path_assembly, rna_id)
            if path_fastq_left == path_fastq_right and path_fastq_left != "-":
                error_nb += 1
                logger.error("l%i: Left and right fastq files are identical, check sample sheet line of %s", l_nb, rna_id)
            for ref in ref_species.split(","):
                if not ref in All_Species:
                    error_nb += 1
                    logger.error("l%i: %s is not in the species tree.\nSpecies in the species tree:\n\t%s", l_nb,  ref, "\n\t".join(All_Species))
            if sp in All_Species:
                if run_apytram.strip() in ["y", "yes", "Y", "Yes"]:
                    RefSpApytram.extend(ref_species.split(","))
                if run_trinity.strip() in ["y", "yes", "Y", "Yes"]:
                    RefSpTrinity.extend(ref_species.split(","))
                if sp not in RnaSp:
                    RnaSp.append(sp)
                if orientation not in ["FR", "RF", "F", "R", "US", "UP","-"]:
                    error_nb += 1
                    logger.error("l%i: orientation must be  in [FR,RF,F,R,US,UP,-] and not: %s", l_nb, orientation)
            else:
                error_nb += 1
                logger.error("l%i: %s is not in the species tree.\nSpecies in the species tree:\n\t%s", l_nb, sp, "\n\t".join(All_Species))

        else:
            logger.error("l%i: Config file has not 11 elements in line:\nid	species	apytram_group	ref_species	path_fastq_single	path_fastq_left	path_fastq_right	orientation	run_trinity	path_assembly	run_apytram)\nbut:\n%s", l_nb, line)
            error_nb += 1


    logger.debug("Ref Trinity : %s", RefSpTrinity)
    logger.debug("Ref Apytram : %s", RefSpApytram)
    logger.debug("Rna species : %s", RnaSp)

if error_nb > 0:
    logger.error("%s errors in the sample sheet. Please correct them.\n", error_nb)
    sys.exit(1)

### Read all files in seq2sp_dir
Seq2Sp_dict = {}

def read_seq2species_file(Seq2Sp_dict, File):
    if os.path.isfile(File):
        f = open(File, "r")
        for line in f:
            if line.strip() == "":
                pass
            else:
                (seq, sp) = line.split("\t")
                if not Seq2Sp_dict.has_key(seq):
                    Seq2Sp_dict[seq] = sp.replace("\n", "")
                else:
                    logger.error("ERROR : Sequence name \"%s\" is not unique", seq)
                    sys.exit(1)
        f.close()
    return Seq2Sp_dict

logger.info("Parse each Sequence-Species link file")
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
                AliDict.setdefault(name, []).append(line.replace("\n", ""))
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
                           '\n'.join(seq[i:i+60] for i in range(0, len(seq), 60)),"\n"])

        with open(Transcriptome_File, "a") as f:
            f.write("".join(string))

        #Tab Seq 2 Fam:
        FamSeqLink_File = "%s/%s_Fam_Seq.tsv" %(SeqFamLinkDirPath, sp)
        string = []

        for name in Ref_dic_trinity[sp]:
            string.extend([name, "\t", Family, "\n"])

        with open(FamSeqLink_File, "a") as f:
            f.write("".join(string))

def write_seq_ref_apytram(Ref_dic_trinity, AliDict_i, Family):
    for sp in Ref_dic_apytram.keys():
        #gene family:
        GeneFamily_File = "%s/%s.%s.fa" %(ApytramGeneFamDirPath, sp, Family)
        string = []
        for name in Ref_dic_trinity[sp]:
            seq = ''.join(AliDict_i[name]).replace("-", "")
            string.extend([">", name, "\n",
                           '\n'.join(seq[i:i+60] for i in range(0, len(seq), 60)),"\n"])

        with open(GeneFamily_File, "w") as f:
            f.write("".join(string))

SeenSeq2SpDict = {}
CountDict = {}
CountDict2 = {}
for sp in All_Species:
    CountDict2.setdefault(sp, {"Nb_seq":0, "Nb_family":0, "Families":[]})
Nb_Family = 0

UsableFam_list = []
NotUsableFam_dict = {}
ErrorFam_dict = {}
Nb_NoError_inFam = 0


logger.info("Parse each fasta file")
for f in glob.glob("%s/*" %ali_dir):
    Family = os.path.basename(f).split('.')[0]
    Extention = os.path.basename(f).split('.')[1]
    print Family
    
    CountDict[Family] = {"Nb_seq": 0, "Nb_sp": 0, "Status": "Usable"}
    Nb_Family += 1
    UsableFam = True
    NoError_inFam = True

    SeenSeq2SpDict_i = {}
    SpeciesList = []
    
    AliDict_i, err = read_ali_file(f)
    Nb_seqs = len(AliDict_i.keys())
    CountDict[Family]["Nb_seq"] = Nb_seqs

    if NoError_inFam and err:
        Reason = "%s is not a fasta file" %(f)
        logger.error("[%s] -->\t%s",Family,Reason)
        ErrorFam_dict.setdefault(Family,[]).append(Reason)
        NoError_inFam = False

    if NoError_inFam:
        for seq in AliDict_i.keys():
            print seq
            if NoError_inFam and seq.replace(sep, "") != seq :
                Reason = "The sequence name [%s] contained the separator [%s] used in CAARS. Please remove it." %(seq, sep)
                ErrorFam_dict.setdefault(Family,[]).append(Reason)
                NoError_inFam = False
            
            #check if the seq is in the sp2seq dir
            if NoError_inFam:
                if seq in Seq2Sp_dict:
                    sp = Seq2Sp_dict[seq]
                    SpeciesList.append(sp)
                else:
                    Reason = "No sequence called %s in the seq2sp-dir" %(seq)
                    ErrorFam_dict.setdefault(Family,[]).append(Reason)
                    NoError_inFam = False
            
            # Check if the sequence is unique
            if NoError_inFam and SeenSeq2SpDict.has_key(seq):
                Reason = "Sequence names are not unique (ex:  %s)" %(seq)
                ErrorFam_dict.setdefault(Family,[]).append(Reason)
                NoError_inFam = False

            SeenSeq2SpDict[seq] = sp
            SeenSeq2SpDict_i[seq] = sp
            
            # Check if the sp is in the species tree
            if NoError_inFam and not sp in All_Species:
                Reason = "%s not in the species tree (%s)" %(sp, species_tree_file)
                ErrorFam_dict.setdefault(Family,[]).append(Reason)
                NoError_inFam = False

    Nb_sp = len(set(SpeciesList))
    CountDict[Family]["Nb_sp"] = Nb_sp

    if  NoError_inFam and Extention != "fa":
        Reason = "%s is not a fasta file with Family.fa as filename. (Detected extention %s)" %(f, Extention)
        ErrorFam_dict.setdefault(Family,[]).append(Reason)
        NoError_inFam = False

    if NoError_inFam and not os.path.isfile("%s/%s.%s" %(ali_dir, Family, "fa")):
        Reason = "%s is not a fasta file with Family.fa as filename.(Detected file: %s/%s.%s)" %(f, ali_dir, Family, "fa")
        ErrorFam_dict.setdefault(Family,[]).append(Reason)
        NoError_inFam = False

    # Check all sequence name in Seq2SpDict
    if NoError_inFam and not len(AliDict_i.keys()) == len(set(AliDict_i.keys()).intersection(set(Seq2Sp_dict.keys()))):
        Reason = "All sequences present in %s are not in a file from %s" %(f, seq2sp_dir)
        ErrorFam_dict.setdefault(Family,[]).append(Reason)
        NoError_inFam = False
    
    if NoError_inFam:
        # Check only ACTG- character in ali and same length
        invalid_char = set()
        invalid_char_seq = []
        length_seq = []
        for n, s in AliDict_i.items():
            seq = "".join(s)
            length_seq.append(len(seq))
            invalid_char_tmp = set(re.sub("[atgcnuwsmkrybdhvATGCNUWSMKRYBDHV-]","", seq))
            if invalid_char_tmp:
                invalid_char |= invalid_char_tmp
                invalid_char_seq.append(n)

        if invalid_char:
            Reason = "Invalid character [%s] present in [%s]." %(",".join(list(invalid_char)), ",".join(invalid_char_seq))
            ErrorFam_dict.setdefault(Family,[]).append(Reason)
            NoError_inFam = False

        length_seq = list(set(length_seq))
        if NoError_inFam and len(length_seq) !=1:
            Reason = "Sequences must be aligned. Different sequence lengths [%s]." %(",".join(map(str,length_seq)))
            ErrorFam_dict.setdefault(Family,[]).append(Reason)
            NoError_inFam = False

    if NoError_inFam and Nb_seqs < 3:
        Reason = "%s has less than 3 sequences. (%s sequences detected in: %s)" %(Family, Nb_seqs, f)
        if Nb_sp < 3:
            Reason += "AND %s has less than 3 species. (%s species detected in: %s)" %(Family, Nb_sp, f)
        logger.warning("[%s] : %s",Family,Reason)
        NotUsableFam_dict.setdefault(Family,[]).append(Reason)
        UsableFam = False

    if NoError_inFam and UsableFam and Nb_sp < 3:
        Reason = "%s has less than 3 species. (%s species detected in: %s)" %(Family, Nb_sp, f)
        logger.warning("[%s] : %s",Family,Reason)
        NotUsableFam_dict.setdefault(Family,[]).append(Reason)
        UsableFam = False

    if NoError_inFam:
        # Write each intermediary files
        Ref_dic_trinity = {}
        Ref_dic_apytram = dict([(key, []) for key in RefSpApytram])

        for seq in AliDict_i.keys():
            sp = Seq2Sp_dict[seq]
            if sp in RefSpTrinity + RefSpApytram:
                #if sp in RefSpTrinity:
                Ref_dic_trinity.setdefault(sp, []).append(seq)
                if sp in RefSpApytram:
                    Ref_dic_apytram.setdefault(sp, []).append(seq)
            # Some log
            CountDict2[sp]["Nb_seq"] += 1
            if not Family in CountDict2[sp]["Families"]:
                CountDict2[sp]["Nb_family"] += 1
                CountDict2[sp]["Families"].append(Family)

    if NoError_inFam:
        write_validated_sp2seq(SeenSeq2SpDict_i, Family)
        write_seq_ref_Trinity(Ref_dic_trinity, AliDict_i, Family)
        write_seq_ref_apytram(Ref_dic_apytram, AliDict_i, Family)
    else:
        Nb_NoError_inFam+=1

    if NoError_inFam and UsableFam:
        UsableFam_list.append(Family)
    else:
        CountDict[Family]["Status"] = "NotUsable"

if Nb_NoError_inFam != 0:
    logger.error("%s families with errors:", Nb_NoError_inFam)
    logger.error("Correct or remove them (See below):")
    for (Family, Reasons) in ErrorFam_dict.items():
        logger.error("[%s] : %s", Family, " ".join(Reasons)) 
    sys.exit(1)

with open(out_dir+"/UsableFamilies.txt", "w") as f:
    f.write("\n".join(sorted(UsableFam_list)) + "\n")

# Statistics:
with open(out_dir+"/SpeciesMetadata.txt", "w") as f:
    f.write("\t".join(["Species","Nb_seq", "Nb_family"])+"\n")
    for sp in All_Species:
        f.write("\t".join(map(str,[sp,CountDict2[sp]["Nb_seq"], CountDict2[sp]["Nb_family"]]))+"\n")

# Statistics:
with open(out_dir+"/FamilyMetadata.txt", "w") as f:
    f.write("\t".join(["Family","Nb_seq", "Nb_sp", "Status"])+"\n")
    for Family in sorted(CountDict.keys()):
        f.write("\t".join(map(str,[Family,CountDict[Family]["Nb_seq"], CountDict[Family]["Nb_sp"], CountDict[Family]["Status"]]))+"\n")

for sp in RefSpTrinity + RefSpApytram:
    if not CountDict2[sp]:
        logger.error("No sequence for  %s", sp)
        sys.exit(1)

sys.exit(0)
