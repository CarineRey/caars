#!/usr/bin/python
# coding: utf-8

# File: PhyldogPrepData.py
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
import argparse

#recursively builds a list "liste" of all the files in a directory "nomdir"
def explore(nomdir, liste):
    listefiles = os.listdir(nomdir)
    for i in listefiles:
        if os.path.isdir(os.path.join(nomdir, i)):
            explore(os.path.join(nomdir, i), liste)
        elif not i.startswith("."):
            liste.append(os.path.abspath(os.path.join(nomdir, i)))
    return()


### Option defining
parser = argparse.ArgumentParser(prog="PhyldogPrepData.py",
                                 description='''
    A script to prepare to lauch phyldog (derived from http://www.prabi.fr/redmine/attachments/download/561/prepareData.py)''')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')


##############
requiredOptions = parser.add_argument_group('Options')
requiredOptions.add_argument('-seq', type=str,
                             help="Absolute PATH to the alignment file", required=True)
requiredOptions.add_argument('-datatype', type=str, choices=["DNA", "RNA", "CODON", "PROTEIN"], default="DNA",
                             help="What alphabet? DNA or RNA or CODON or PROTEIN?")
requiredOptions.add_argument('-dataformat', type=str, choices=["FASTA", "MASE", "PHYLIP", "CLUSTAL", "DCSE", "NEXUS"], default="FASTA",
                             help="What alignment format? FASTA or MASE or PHYLIP or CLUSTAL or DCSE or NEXUS?")

requiredOptions.add_argument('-starting_tree', type=str,
                             help="Absolute PATH to the starting tree file")

requiredOptions.add_argument('-link', type=str,
                             help="Absolute PATH to the link file", required=True)

requiredOptions.add_argument('-optdir', type=str,
                             help="Absolute PATH to the folder that will contain all the option files", required=True)

requiredOptions.add_argument('-species_tree_resdir', type=str, default="",
                             help="Absolute PATH to the folder that will contain result files concerning the species tree")
requiredOptions.add_argument('-gene_trees_resdir', type=str, default="",
                             help="Absolute PATH to the folder that will contain result files concerning all gene trees")

requiredOptions.add_argument('-species_tree_file', type=str,
                             help="Absolute PATH to the starting species tree file (Newick)")
requiredOptions.add_argument('-startingtree', type=str, choices=["RANDOM", "MRP"], default="RANDOM",
                             help="If no species tree file, What kind of species tree? RANDOM or MRP?")

requiredOptions.add_argument('-topospecies', action='store_true', default=False,
                             help="If you want to optimize the species tree topology")

requiredOptions.add_argument('-dlopt', type=str, choices=["BRANCHWISE", "AVERAGE"], default="",
                             help="If you want to optimize the duplication and loss parameters? Branchwise or average ?")


requiredOptions.add_argument("-equgenomes", action='store_true', default=False,
                             help="Do you want to assume that all genomes, if they had no missing data, would have the same number of genes?")
requiredOptions.add_argument("-topogene", action='store_true', default=False,
                             help="Do you want to optimize the gene tree topologies?")

requiredOptions.add_argument('-timelimit', type=int, default=2,
                             help="What time limit do you want to set, in hours (minimum 2h)?")


##############

### Option parsing
args = parser.parse_args()

print " ".join(sys.argv)

SEQ = args.seq
DATATYPE = args.datatype
DATAFORMAT = args.dataformat
LINK = args.link
OPTDIR = args.optdir
SP_RESDIR = args.species_tree_resdir
GENES_RESDIR = args.gene_trees_resdir
STARTINGTREE = args.starting_tree
TREEFILE = args.species_tree_file.replace("//","/")
STARTINGTREE = args.startingtree
if TREEFILE:
    TREEFILEGIVEN = True
else:
    TREEFILEGIVEN = False

TOPOSPECIES = args.topospecies

DLOPT = args.dlopt
if args.dlopt:
    DLPARAM = True
else:
    DLPARAM = False

EQUGENOMES = args.equgenomes
TOPOGENE = args.topogene
TIMELIMIT = str(args.timelimit)

if __name__ == '__main__':
    #Alignment

    if not os.path.isfile(SEQ):
        print SEQ + " does not exist.\n"
        sys.exit()
    listAlns = [os.path.abspath(SEQ)]
    #StartingTree

    if STARTINGTREE:
        if not os.path.isfile(STARTINGTREE):
            print STARTINGTREE + " does not exist.\n"
            sys.exit()
    listTrees = [os.path.abspath(STARTINGTREE)]
    #Data type
    if not DATATYPE in ["DNA", "RNA", "CODON", "PROTEIN"]:
        print DATATYPE + " is not an option I understand for an alphabet.\n"
        sys.exit()

    #File format
    if not DATAFORMAT in ["FASTA", "MASE", "PHYLIP", "CLUSTAL", "DCSE", "NEXUS"]:
        print DATAFORMAT + " is not an option I understand for an alignment format.\n"
        sys.exit()

    #link file folder

    pathExists = os.path.isfile(LINK)
    if not pathExists:
        print LINK + " does not exist."
        sys.exit()

    listLinks = [os.path.abspath(LINK)]

    #Option file folder
    pathExists = os.path.isdir(OPTDIR)
    if not pathExists:
        print "\nCreating the folder " + OPTDIR + "\n"
        try:
            os.makedirs(OPTDIR)
        except OSError:
            if not os.path.isdir(OPTDIR):
                raise
    print OPTDIR + " exists and will contain the option files.\n"

    #Result file folder
    for RESDIR in [SP_RESDIR, GENES_RESDIR]:
        if RESDIR:
            if not os.path.isdir(RESDIR):
                print "\nCreating the folder "+RESDIR+"\n"
                try:
                    os.makedirs(RESDIR)
                except OSError:
                    if not os.path.isdir(RESDIR):
                        raise
            print RESDIR + " exists and will contain the result files.\n"

    #Species tree
    if TREEFILEGIVEN:
        pathExists = os.path.isfile(TREEFILE)
        if not pathExists:
            sys.exit(1)
        print "Starting from "+TREEFILE+"\n"
    else:
        if not STARTINGTREE in ["RANDOM", "MRP"]:
            sys.exit(1)
        if STARTINGTREE == "RANDOM":
            print "Starting species tree: RANDOM\n"
        else:
            print "Starting species tree: MRP\n"


    #Do we optimize the species tree?
    if TOPOSPECIES:
        print "Optimizing the species tree topology.\n"
    else:
        print "NOT optimizing the species tree topology.\n"

    #Do we optimize the duplication and loss parameters, and how?
    if DLPARAM:
        if not DLOPT in ["BRANCHWISE", "AVERAGE"]:
            if DLOPT == "BRANCHWISE":
                print "Optimizing the duplication and loss parameters, using BRANCHWISE.\n"
            else:
                print "Optimizing the duplication and loss parameters, using AVERAGE.\n"
    else:
        print "NOT optimizing the duplication and loss parameters.\n"

    #Do we want to assume that all species have roughly the same number of genes,
    #and use that for computing an expected amount of missing data?
    if EQUGENOMES:
        print "Assuming that all the genomes have the same number of genes, so that all missing genes are artifactual.\n"
    else:
        print "NOT assuming that all the genomes have the same number of genes.\n"

    #Do we optimize the gene trees?
    if TOPOGENE:
        print "Optimizing the gene tree topologies.\n"
    else:
        print "NOT optimizing the gene tree topologies (but they will still be rerooted).\n"

    #Time limit

    print "We will stop at the latest about 1 hour before %s.\n" %TIMELIMIT

    ###########################################
    ###########################################
    #Creating gene family-specific option files
    listSpecies = list()
    listOptionFiles = list()
    listSizes = list()
    dictLinks = dict()
    dictTrees = dict()
    print "Now we create one option file per gene family, but we do so only if we find an alignment file and a link file that correspond to each other.\n"
    print "Correspondance is based on the radical of the file names (file name without .extension).\n"
    listAlns.sort()
    listLinks.sort()
    listTrees.sort()
    radical = ""
    for link in listLinks:
        dictLinks[os.path.basename(link).split('.')[0]] = link
    if STARTINGTREE:
        for tree in listTrees:
            dictTrees[os.path.basename(tree).split('.')[0]] = tree

    for aln in listAlns:
        radical = os.path.basename(aln).split('.')[0]
        if dictLinks.__contains__(radical):
            #First, open the link file in order to build a list of species
            try:
                fin = open(dictLinks[radical], 'r')
            except IOError, e:
                print "Unknown file: ", dictLinks[radical]
                sys.exit()
            for l in fin:
                num = 1
                sp = l.split(":")[0]
                if ", " in l:
                    num = len(l.split(";"))
                for i in range(num):
                    listSpecies.append(sp)
            fin.close()
            try:
                fopt = open(os.path.join(OPTDIR, radical)+'.opt', 'w')
            except IOError, e:
                print "Unknown file: ", os.path.join(OPTDIR, radical) + '.opt'
                sys.exit()
            fopt.write("\n######## First, data files ########\n")
            #fopt.write("PATH="+os.path.join(SEQDIR, "")+"\n")
            if GENES_RESDIR:
                fopt.write("RESULT="+os.path.join(GENES_RESDIR, "")+"\n")
            else:
                fopt.write("RESULT=\n")
            fopt.write("DATA="+radical+"\n")
            fopt.write("taxaseq.file="+dictLinks[radical]+"\n")
            fopt.write("input.sequence.file="+aln+"\n")
            if "FASTA".startswith(DATAFORMAT):
                fopt.write("input.sequence.format=Fasta\n")
            elif "MASE".startswith(DATAFORMAT):
                fopt.write("input.sequence.format=Fasta\n")
            elif "PHYLIP".startswith(DATAFORMAT):
                fopt.write("input.sequence.format=Phylip\n")
            elif "CLUSTAL".startswith(DATAFORMAT):
                fopt.write("input.sequence.format=Clustal\n")
            elif "DCSE".startswith(DATAFORMAT):
                fopt.write("input.sequence.format=Dcse\n")
            elif "NEXUS".startswith(DATAFORMAT):
                fopt.write("input.sequence.format=Nexus\n")
            listSizes.append(str(os.stat(aln)[6]))
            if dictTrees.__contains__(radical):
                fopt.write("gene.tree.file=" + dictTrees[radical]+ "\n")
                fopt.write("init.gene.tree=user\n")
            else:
                fopt.write("gene.tree.file=$(RESULT)$(DATA).GeneTree\n")
                fopt.write("init.gene.tree=bionj\n") #the program will build a starting user gene tree
                if STARTINGTREE:
                    print "Alignment "+ aln + " does not have a corresponding starting tree file.\n"
            fopt.write("output.reconciled.tree.file=$(RESULT)$(DATA).ReconciledTree\n")
            fopt.write("output.duplications.tree.file=$(RESULT)$(DATA).DuplicationTree\n")
            fopt.write("output.losses.tree.file=$(RESULT)$(DATA).LossTree\n")
            fopt.write("output.numbered.tree.file=$(RESULT)$(DATA).NumberedTree\n")
            fopt.write("input.sequence.sites_to_use=all\n")
            fopt.write("input.sequence.max_gap_allowed=100%\n")
            fopt.write("output.starting.gene.tree.file=$(RESULT)$(DATA).StartingTree\n")
            fopt.write("\n######## Second, model options ########\n")
            if DATATYPE == "DNA":
                fopt.write("alphabet=DNA\n")
                fopt.write("model=GTR(a=1.17322, b=0.27717, c=0.279888, d=0.41831, e=0.344783, theta=0.523374, theta1=0.542411, theta2=0.499195)\n")
                fopt.write("rate_distribution=Gamma(n=4,alpha=1)\n")
            elif DATATYPE == "RNA":
                fopt.write("alphabet=RNA\n")
                fopt.write("model=GTR( initFreqs=observed )\n")
                fopt.write("rate_distribution=Gamma(n=4,alpha=1)\n")
            elif DATATYPE == "CODON":
                fopt.write("alphabet=Codon(letter=DNA)\n")
                fopt.write("input.sequence.remove_stop_codons = yes\n")
                fopt.write("genetic_code=Standard\n")
                fopt.write("model=YN98( kappa=1, omega=1.0, initFreqs=observed )\n")
            elif DATATYPE == "PROTEIN":
                fopt.write("alphabet=Protein\n")
                fopt.write("model=LG08+F(initFreqs=observed )\n")
                fopt.write("rate_distribution=Gamma(n=4,alpha=1)\n")
            #fopt.write("optimization.ignore_parameter=dist_Gamma.alpha, GTR.a, GTR.b, GTR.c, GTR.d, GTR.e, GTR.theta, GTR.theta1, GTR.theta2\n")
            #fopt.write("\n######## Then, algorithm options ########\n")
            #fopt.write("heuristics.level=0\n")
            #fopt.write("species.id.limit.for.root.position=3 #Useless unless heuristics.level=1\n")
            fopt.write("\n######## Finally, optimization options ########\n")
            if TOPOGENE:
                fopt.write("optimization.topology=yes\n")
            else:
                fopt.write("optimization.topology=no\n")
            #fopt.write("optimization.topology.algorithm_nni.method=fast\n")
            #fopt.write("optimization.tolerance=0.01\n")
            #fopt.write("optimization.method_DB.nstep=0\n")
            #fopt.write("optimization.topology.numfirst=false\n")
            #fopt.write("optimization.topology.tolerance.before=100\n")
            #fopt.write("optimization.topology.tolerance.during=100\n")
            #fopt.write("optimization.max_number_f_eval=1000000\n")
            #fopt.write("optimization.final=none\n")
            #fopt.write("optimization.verbose=0\n")
            #fopt.write("optimization.message_handler=none\n")
            #fopt.write("optimization.profiler=none\n")
            #fopt.write("optimization.reparametrization=no\n")
            fopt.close()
            listOptionFiles.append(os.path.join(OPTDIR, radical) + '.opt')
        else:
            print "Alignment " + aln + " does not have a corresponding link file. Skipping it.\n"



    #Now the gene option files have been created

    #Creating the list of gene option files
    try:
        fout = open(os.path.join(OPTDIR, "listGenes.txt"), 'w')
    except IOError, e:
        print "Unknown file: ", os.path.join(OPTDIR, "listGenes.txt")
        sys.exit()
    for i in range(len(listOptionFiles)):
        fout.write(listOptionFiles[i] + ":" + listSizes[i] + "\n")
    fout.close()

    #Creating the species list file.
    dictSpecies = dict()
    for i in listSpecies:
        if dictSpecies.__contains__(i):
            dictSpecies[i] = dictSpecies[i] + 1
        else:
            dictSpecies[i] = 1

    #We have some data, why not output it?
    print "\n\t~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
    print "\n\tTaxonomic distribution in the gene families for which we created option files:\n"
    max = 0
    min = 1000000
    sum = 0
    try:
        fsp = open(os.path.join(OPTDIR, "listSpecies.txt"), 'w')
    except IOError, e:
        print "Unknown file: ", os.path.join(OPTDIR, "listSpecies.txt")
        sys.exit()
    for (k, v) in dictSpecies.items():
        print k + " : " + str(v) + " genes."
        fsp.write(k+"\n")
        sum = sum +v
        if v > max:
            max = v
        if v < min:
            min = v
    fsp.close()
    if len(dictSpecies) > 0:
        print "\n\tAverage number of genes: "+str(sum/len(dictSpecies))+" ; Maximum: "+str(max)+" ; Minimum: "+str(min)+"\n"
    else:
        print "\n\tNo gene\n"


    #Creating a taxonomic coverage file
    if EQUGENOMES:
        try:
            fsp = open(os.path.join(OPTDIR, "listSpeciesWithSequenceCoverage.txt"), 'w')
        except IOError, e:
            print "Unknown file: ", os.path.join(OPTDIR, "listSpeciesWithSequenceCoverage.txt")
            sys.exit()
        species = sorted(dictSpecies.keys())
        for sp in dictSpecies.keys():
            v = dictSpecies[sp]
            fsp.write(sp + " : " + str(100*v/max) + " \n")
        fsp.close()


    #Now, the general options.
    try:
        fopt = open(os.path.join(OPTDIR, "GeneralOptions.txt"), 'w')
    except IOError, e:
        print "Unknown file: ", os.path.join(OPTDIR, "GeneralOptions.txt")
        sys.exit()
    fopt.write("\n######## First, data files ########\n")
    fopt.write("OPT="+os.path.join(OPTDIR, "")+"\n")
    if SP_RESDIR:
        fopt.write("RESULT="+os.path.join(SP_RESDIR, "")+"\n")
    else:
        fopt.write("RESULT=\n")
    fopt.write("DATA="+radical+"\n")
    if TREEFILEGIVEN:
        fopt.write("init.species.tree=user\n")
        fopt.write("species.tree.file="+TREEFILE+"\n")
    elif STARTINGTREE == "RANDOM":
        fopt.write("init.species.tree=random #user\n")
    else:
        fopt.write("init.species.tree=mrp\n")
    fopt.write("species.names.file=$(OPT)listSpecies.txt\n")
    fopt.write("starting.tree.file=$(RESULT)StartingTree.tree\n")
    fopt.write("output.tree.file=$(RESULT)OutputSpeciesTree.tree\n")
    #There follows a list of genes for which a tree needs to be built.
    fopt.write("genelist.file=$(OPT)listGenes.txt\n")
    fopt.write("output.duplications.tree.file=$(RESULT)OutputSpeciesTree_ConsensusDuplications.tree\n")
    fopt.write("output.losses.tree.file=$(RESULT)OutputSpeciesTree_ConsensusLosses.tree\n")
    fopt.write("output.numbered.tree.file=$(RESULT)OutputSpeciesTree_ConsensusNumbered.tree\n")
    fopt.write("\n######## Second, options ########\n")
    if TOPOSPECIES:
        fopt.write("optimization.topology=yes\n")
    else:
        fopt.write("optimization.topology=no\n")
    #fopt.write("species.id.limit.for.root.position=3\n")
    if DLPARAM:
        if DLOPT == "BRANCHWISE":
            fopt.write("branchProbabilities.optimization=average_then_branchwise\n")
        else:
            fopt.write("branchProbabilities.optimization=average\n")
    else:
        fopt.write("branchProbabilities.optimization=no")
    if EQUGENOMES:
        fopt.write("genome.coverage.file=$(PATH)HomolensSpeciesSequenceCoverage\n")
    fopt.write("spr.limit=5\n")
    fopt.write("time.limit="+TIMELIMIT+"\n")
    fopt.close()

    print "\n\t~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
    print "\n\n\n\tTo launch PHYLDOG, you will need to run something like:\n\n"

    print "\tmpirun -np NUMBER_OF_PROCESSES PATH_TO_PHYLDOG/phyldog param=%s\n\n" %(os.path.join(OPTDIR, "GeneralOptions.txt"))

    print "\n\tThank you for using this PHYLDOG file preparation script!\n\n"



