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
import subprocess
import logging

class Fasttree(object):
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


    def get_output(self, output=""):
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
            command.extend(["-out", output])
        elif self.OutputTree:
            command.extend(["-out", self.OutputTree])

        command.append(self.InputAliFile)
        self.logger.debug(" ".join(command))
        try:
            Out = subprocess.call(command,
                                  stdout=open("/dev/null", "w"),
                                  stderr=open("/dev/null", "w"))
        except:
            os.system(
            "echo Unexpected error when we launch fasttree:\n"
            )
            print " ".join(command)
        return Out

class Phylomerge(object):
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


    def launch(self, output=""):
        command = ["phylomerge",
                   "input.sequence.file=%s" %self.InputAliFile,
                   "input.method=%s" %self.InputMethod,
                   "choice.criterion=%s" %self.ChoiceCriterion,
                   "deletion.method=%s" %self.DeletionMethod
                  ]

        if output:
            self.OutputSequenceFile = output

        if self.OutputSequenceFile:
            command.append(
            "output.sequence.file=%s" %self.OutputSequenceFile
            )

        if self.InputTreeFile:
            if os.path.isfile(self.InputTreeFile):
                command.append(
                "input.tree.file=%s" %self.InputTreeFile
                )

        if self.TaxonsToRefine:
            if os.path.isfile(self.TaxonsToRefine):
                command.append(
                "taxons.to.refine=%s" %self.TaxonsToRefine
                )

        if self.TaxonToSequence:
            if os.path.isfile(self.TaxonToSequence):
                command.append(
                "taxon.to.sequence=%s" %self.TaxonToSequence
                )
        elif self.SequenceToTaxon:
            if os.path.isfile(self.SequenceToTaxon):
                command.append(
                "sequence.to.taxon=%s" %self.SequenceToTaxon
                )

        if self.PrescreeningOnSizeByTaxon:
            command.append("prescreening.on.size.by.taxon=yes")
        else:
            command.append("prescreening.on.size.by.taxon=no")

        if self.SelectionByTaxon:
            command.append("selection.by.taxon=yes")

        if self.RearrangeTree:
            command.append("rearrange.tree=yes")

        if self.BootstrapThreshold:
            command.append(
            "bootstrap.threshold=%s" %self.BootstrapThreshold
            )

        if self.OutputTaxonToSequence:
            command.append(
            "output.taxon.to.sequence=%s" %self.OutputTaxonToSequence
            )

        self.logger.debug(" ".join(command))

        p = subprocess.Popen(command,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        (out, err) = p.communicate()
        if err:
            self.logger.error(err)
        return (out, err)
