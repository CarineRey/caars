#!/usr/bin/python
# coding: utf-8

# File: annote_reconciledtree.py
# Created by: Carine Rey
# Created on: march 2010
#
#
# Copyright 2017 Carine Rey
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
import re
import sys
import time
import shutil
import logging
import argparse

from ete3 import Tree

start_time = time.time()


### Option defining
parser = argparse.ArgumentParser(prog="annote_reconciledtree.py",
                                 description='''Detect D and L event in a reconciled tree''')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')


##############
requiredOptions = parser.add_argument_group('Required arguments')
requiredOptions.add_argument('-s', '--sp_tree', type=str,
                             help='Species tree filename.', required=True)
requiredOptions.add_argument('-t', '--tree', type=str,
                             help='Gene tree filename.', required=True)
requiredOptions.add_argument('-sp2seq', type=str,
                             help='Link file name. (spA:seq1,seq2,...', required=True)
requiredOptions.add_argument('-out', '--output_prefix', type=str, default="./output",
                   help="Output prefix (Default ./output)")
##############


##############
Options = parser.add_argument_group('Options')
Options.add_argument('-log', type=str, default="annote_reconciledtree.log",
                   help="a log file to report avancement (default: annote_reconciledtree.log)")
Options.add_argument('--debug', action='store_true', default=False,
                   help="debug mode, default False")

### Option parsing
args = parser.parse_args()

### Read arguments
SpeciesTreeFilename = args.sp_tree
GeneTreeFilename    = args.tree
Sp2SeqFilename      = args.sp2seq

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
def end(ReturnCode):
    logger.debug("--- %s seconds ---", str(time.time() - start_time))
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

for inputfile in [SpeciesTreeFilename, GeneTreeFilename, Sp2SeqFilename]:
    if not os.path.isfile(inputfile):
        logger.error(inputfile +" is not a file.")
        end(1)


from ete3 import PhyloTree

#read sp2seq file and build a map dictionnary


with open(Sp2SeqFilename,"r") as File:
    sp2seq_list = File.read().strip().split("\n")

seq2sp_dict = {}
for sp_seqs in sp2seq_list:
    sp, seqs = sp_seqs.split(":")
    for seq in seqs.split(","):
        seq2sp_dict[seq]=sp


def get_species_name(node_name_string):
    return seq2sp_dict[node_name_string]
    
def put_species_name(node_name_string):
    return node_name_string
    

# read the gene tree
genetree = PhyloTree(GeneTreeFilename, sp_naming_function=get_species_name)
sptree = PhyloTree(SpeciesTreeFilename, sp_naming_function=put_species_name)

logger.debug("Genetree")
for n in genetree.get_leaves():
    logger.debug("node: %s Species name: %s", n.name, n.species)

logger.debug("SpeciesTree")
for n in sptree.get_leaves():
    logger.debug("node: %s Species name: %s", n.name, n.species)

iS = 0
sp_dict = {}
for n in sptree.traverse():
    n.S=iS
    iS+=1
    if not n.is_leaf():
        n.name = n.S
    else:
        sp_dict[n.name] = n.S


logger.debug("sp_dict")
logger.debug(sp_dict)

# Let's reconcile our genetree with the species tree
recon_tree, events = genetree.reconcile(sptree)
# a new "reconcilied tree" is returned. As well as the list of
# inferred events.


ntrees, ndups, sptrees =  genetree.get_speciation_trees()
logger.debug( "Found %d species trees and %d duplication nodes", ntrees, ndups)

HomologySummary=["DUPLICATIONS : " + str(ndups) ]

logger.debug( "Orthology and Paralogy relationships:")
for ev in events:
    if ev.etype == "S":
        logger.debug("".join(['ORTHOLOGY RELATIONSHIP:', ','.join(ev.inparalogs), " <===> ", ','.join(ev.orthologs)]))
        HomologySummary.append("".join(['ORTHOLOGY RELATIONSHIP: ', ', '.join(ev.inparalogs), " <===> ", ','.join(ev.orthologs)]))
    elif ev.etype == "D":
        logger.debug("".join(['PARALOGY RELATIONSHIP: ', ', '.join(ev.inparalogs), " <===>" , ','.join(ev.outparalogs)]))
        HomologySummary.append("".join(['PARALOGY RELATIONSHIP: ', ', '.join(ev.inparalogs), " <===> ", ','.join(ev.outparalogs)]))

HomologyFile = OutPrefixName + ".orthologs.txt"
with open(HomologyFile,"w") as File:
        File.write("\n".join(HomologySummary)+"\n")


EventSummary = []
i=0
for n in recon_tree.traverse("postorder"):
    n.ND = i
    if n.is_leaf():
        n.S = sp_dict[n.species]
    if "evoltype" in dir(n):
        n.Ev = n.evoltype
        if n.evoltype == "L":
            EventSummary.append("event(%i,loss)" %(n.S))
        elif n.evoltype == "D":
            sp_dup = n.get_species()
            oldest_sp = sptree.get_common_ancestor(sp_dup)
            n.S = oldest_sp.S
            logger.debug("sp_dup: %s ",sp_dup)
            EventSummary.append("event(%i,duplication)" %(n.S))
    else:
        n.Ev = "S"
    logger.debug("name: %s",n.name)
    logger.debug("S: %s",n.S)
    logger.debug("Ev: %s",n.Ev)
    logger.debug("ND: %s",n.ND)
    i+=1


EventsFile = OutPrefixName + ".events.txt"
with open(EventsFile,"w") as File:
        File.write("\n".join(EventSummary)+"\n")


recon_tree.prune(genetree.get_leaf_names(),preserve_branch_length=True)

i=0
node_2events_and_sp = {}
for n in recon_tree.traverse("postorder"):
    n.ND = i
    node_2events_and_sp[n.ND]={"S": n.S, "Ev":n.Ev}
    i+=1

logger.debug(node_2events_and_sp)


i=0
for n in genetree.traverse("postorder"):
    n.ND=i
    n.S=node_2events_and_sp[n.ND]["S"]
    n.Ev=node_2events_and_sp[n.ND]["Ev"]
    i+=1

genetree.write(format=1, features=[ "Ev", "S", "ND"],outfile=OutPrefixName+"_ReconciledTree.nhx", format_root_node=True)
end(0)
