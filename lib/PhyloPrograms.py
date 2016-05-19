# File: PhyloPrograms.py
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
import subprocess
import logging

class Fasttree:
    """Define an object to launch Fasttree"""
    def __init__(self, InputAliFile):
        self.logger = logging.getLogger('main.lib.Fasttree')
        self.logger.info('creating an instance of Fasttree')
        self.InputAliFile = InputAliFile
        self.OutputTree = ""
        self.QuietOption = False
        self.Gtr = False
        self.Nt = False
        self.Gamma = False
        
        
    def get_output(self, output = ""):
        Out = ""
        command = ["fasttree"]

        if self.Nt:
            command.append("-nt")
        if self.Gtr:
            command.append("-gtr")
        if self.Gamma:
            command.append("-gamma")
        if self.QuietOption:
            command.append("-quiet")     
            
        if output:
            command.extend(["-out",output])
        elif self.OutputTree:
            command.extend(["-out",self.OutputTree])
     
        command.append(self.InputAliFile)
        self.logger.debug(" ".join(command))
        try:
            Out = subprocess.call(command,
                                  stdout=open("/dev/null", "w"),
                                  stderr=open("/dev/null", "w"))
        except:
            os.system("echo Unexpected error when we launch fasttree:\n")
            print " ".join(command)  
        return Out


#Usage for FastTree version 2.1.7 SSE3:
#  FastTree protein_alignment > tree
#  FastTree < protein_alignment > tree
#  FastTree -out tree protein_alignment
#  FastTree -nt nucleotide_alignment > tree
#  FastTree -nt -gtr < nucleotide_alignment > tree
#  FastTree < nucleotide_alignment > tree
#FastTree accepts alignments in fasta or phylip interleaved formats
#
#Common options (must be before the alignment file):
#  -quiet to suppress reporting information
#  -nopr to suppress progress indicator
#  -log logfile -- save intermediate trees, settings, and model details
#  -fastest -- speed up the neighbor joining phase & reduce memory usage
#        (recommended for >50,000 sequences)
#  -n <number> to analyze multiple alignments (phylip format only)
#        (use for global bootstrap, with seqboot and CompareToBootstrap.pl)
#  -nosupport to not compute support values
#  -intree newick_file to set the starting tree(s)
#  -intree1 newick_file to use this starting tree for all the alignments
#        (for faster global bootstrap on huge alignments)
#  -pseudo to use pseudocounts (recommended for highly gapped sequences)
#  -gtr -- generalized time-reversible model (nucleotide alignments only)
#  -wag -- Whelan-And-Goldman 2001 model (amino acid alignments only)
#  -quote -- allow spaces and other restricted characters (but not ' ) in
#           sequence names and quote names in the output tree (fasta input only;
#           FastTree will not be able to read these trees back in)
#  -noml to turn off maximum-likelihood
#  -nome to turn off minimum-evolution NNIs and SPRs
#        (recommended if running additional ML NNIs with -intree)
#  -nome -mllen with -intree to optimize branch lengths for a fixed topology
#  -cat # to specify the number of rate categories of sites (default 20)
#      or -nocat to use constant rates
#  -gamma -- after optimizing the tree under the CAT approximation,
#      rescale the lengths to optimize the Gamma20 likelihood
#  -constraints constraintAlignment to constrain the topology search
#       constraintAlignment should have 1s or 0s to indicates splits
#  -expert -- see more options


class Phylomerge:
    """Define an object to launch Phylomerge"""
    def __init__(self, InputAliFile, InputTreeFile):
        self.logger = logging.getLogger('main.lib.Phylomerge')
        self.logger.info('creating an instance of Phylomerge')
        self.InputAliFile = InputAliFile
        self.InputTreeFile = InputTreeFile
        self.InputMethod = "tree" # "tree" or "matrix"
        self.DeletionMethod = "taxon" # 'threshold' or 'random' or 'sample' or 'taxon'
        self.ChoiceCriterion = "merge" #'length' ou  'length.complete' ou 'merge'
        self.SelectionByTaxon = True
        self.TaxonToSequence = ""
        self.SequenceToTaxon = ""
        self.TaxonsToRefine = ""
        self.TaxonsToRemove = ""
        self.PrescreeningOnSizeByTaxon = False
        self.RearrangeTree = False
        self.BootstrapThreshold = 0
        self.OutputSequenceFile = ""
        self.OutputTaxonToSequence = ""
        
        
    def launch(self, output = ""):
        Out = ""
        command = ["phylomerge", "input.sequence.file=%s" %self.InputAliFile, "input.method=%s" %self.InputMethod,
                   "choice.criterion=%s" %self.ChoiceCriterion, "deletion.method=%s" %self.DeletionMethod
                  ]
        
        if output:
            self.OutputSequenceFile = output
        
        if self.OutputSequenceFile:
            command.append("output.sequence.file=%s" %self.OutputSequenceFile)

        if self.InputTreeFile:
            if os.path.isfile(self.InputTreeFile):
                command.append("input.tree.file=%s" %self.InputTreeFile)
                
        if self.TaxonsToRefine:
            if os.path.isfile(self.TaxonsToRefine):
                command.append("taxons.to.refine=%s" %self.TaxonsToRefine)
                
        if self.TaxonToSequence:
            if os.path.isfile(self.TaxonToSequence):
                command.append("taxon.to.sequence=%s" %self.TaxonToSequence)
        elif self.SequenceToTaxon:
            if os.path.isfile(self.SequenceToTaxon):
                command.append("sequence.to.taxon=%s" %self.SequenceToTaxon)
        
        if self.PrescreeningOnSizeByTaxon:
            command.append("prescreening.on.size.by.taxon=yes")
        else:
            command.append("prescreening.on.size.by.taxon=no")
            
        if self.SelectionByTaxon:
            command.append("selection.by.taxon=yes")
        
        if self.RearrangeTree:
            command.append("rearrange.tree=yes")
        
        if self.BootstrapThreshold:
            command.append("bootstrap.threshold=%s" %self.BootstrapThreshold)
            
        if self.OutputTaxonToSequence:
            command.append("output.taxon.to.sequence=%s" %self.OutputTaxonToSequence)

        self.logger.debug(" ".join(command))
        
        p = subprocess.Popen(command,
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE)
        out, err = p.communicate()
        if err:
            self.logger.error(err)
        return (out,err)


#******************************************************************
#*        Bio++ Phylogenetic Merger/Sampler, version 0.2          *
#* Author: J. Dutheil, B. Boussau            Last Modif. 21/02/14 *
#******************************************************************
#
#__________________________________________________________________________
#phylomerge parameter1_name=parameter1_value
#      parameter2_name=parameter2_value ... param=option_file
#
#      Options considered: 
# - input.method='tree' (could also be a distance matrix with the option 'matrix')
# - input.tree.file=file with the tree in it.
# - deletion.method='threshold' ou 'random' ou 'sample' ou 'taxon'. 'threshold' removes sequences that are so close in the tree that their distance is lower than the 'threshold' value (which is given as another option to the program, default is 0.01). 'sample': random choice of sample_size sequences (default is 10). 'taxon': choice is guided by the identity of the species the sequences come from. In cases several sequences from the same species are monophyletic, a choice will be made according to the 'choice.criterion' option
# - choice.criterion='length' ou  'length.complete' ou 'merge'. 'length' means the longest sequence is selected. 'length.complete' : means the largest number of complete sites (no gaps). 'merge' means that the set of monophyletic sequences is used to build one long 'chimera' sequence corresponding to the merging of them.
# - selection.by.taxon='no' ou 'yes'
# - sequence.to.taxon=linkfile: format: sequence name: species name. Can be replaced by the option taxon.to.sequence
# - taxon.to.sequence=linkfile: format: species name: sequence name. Can be replaced by the option sequence.to.taxon
# - taxons.to.remove= file containing set of species from which sequences should be removed
# - taxons.to.refine= file containing set of species on which the sampling/merging should be done. If not specified, all species are concerned.
# - prescreening.on.size.by.taxon='no' : removes the sequences that are very short compared to other sequences of the same species. If there is only one sequence in this species, it is not removed.
# - output.sequence.file=refined alignment
#__________________________________________________________________________
