#!/usr/bin/python
# coding: utf-8

# File: build_bioindex_fasta.py
# Created by: Carine Rey
# Created on: March 2017
#
#
# Copyright or © or Copr. Carine Rey
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
#


import os
import logging
import time
import sys

from Bio import SeqIO
start_time = time.time()

### Set up the logger
# create logger with 'spam_application'
logger = logging.getLogger('BuildBioIndex')
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

if len(sys.argv) != 3:
    logger.error("2 arguments are required (index, file fasta file")
    sys.exit(1)

index_file = sys.argv[1]
fasta_file = sys.argv[2]


logger.debug("Index file: %s", index_file)
logger.debug("Fasta file: %s", fasta_file)

if os.path.isfile(fasta_file):
    IndexDB = SeqIO.index_db(index_file, fasta_file, "fasta")
else:
    logger.error("Fasta file (%s) is not a file", fasta_file)
    sys.exit(1)




