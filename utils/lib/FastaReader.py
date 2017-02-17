#!/usr/bin/python
# coding: utf-8

# File: FastaReader.py
# Created by: Carine Rey
# Created on: February 2017
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



def write_in_file(String,Filename,mode = "w"):
    if mode in ["w","a"]:
        with open(Filename,mode) as File:
            File.write(String)


class Fasta:
    def __init__(self):
        self.d = {}

    def __str__(self):
        string = []
        for s in self.Sequences:
            string.extend([str(s)])
        return("".join(string))

    def append(self,new_sequence):
        assert isinstance(new_sequence, Sequence), "Sequence must belong to the Sequence class"
        if new_sequence.Complement and new_sequence.Sequence:
            new_sequence.Sequence = complement(new_sequence.Sequence)
            new_sequence.Complement = False
            
        self.d[new_sequence.Name] = new_sequence
    
    def get(self, name, default=""):
        if name in self.d.keys():
            return self.d[name].Sequence
        else:
            return default

    def read_fasta(self, FastaFilename = "" , String = ""):
        if String:
            Fasta = String.strip().split("\n")
        elif os.path.isfile(FastaFilename):
            with open(FastaFilename,"r") as File:
                Fasta = File.read().strip().split("\n")
        else:
            Fasta = []

        name = ""
        sequence = ""
        sequence_list = []

        for line in Fasta + [">"]:
            if re.match(">",line):
                # This is a new sequence write the previous sequence if it exists
                if sequence_list:
                    new_sequence = Sequence()
                    new_sequence.Name = name
                    new_sequence.Sequence = "".join(sequence_list)
                    self.append(new_sequence)
                    sequence_list = []

                name = line.split()[0][1:] # remove the >

            elif name != "":
                sequence_list.append(line)
            else:
                pass

    def filter_fasta(self, SelectedNames):
        FilteredFasta = Fasta()
        for s in self.Sequences:
            if s.Name in SelectedNames:
                FilteredFasta.append(s)
        return FilteredFasta

    def dealign_fasta(self):
        DealignedFasta = Fasta()
        for s in self.Sequences:
            s.Sequence = s.Sequence.replace("-", "")
            DealignedFasta.append(s)
        return DealignedFasta

    def write_fasta(self, OutFastaFile):
        # Write all sequences in the file
        write_in_file(str(self), OutFastaFile)


class Sequence:
    def __init__(self,):
        self.Name = ""
        self.Sequence = ""
        
        ### Caracteristics
        self.Complement = False

    def __str__(self):
        return(">" + self.Name + "\n" + '\n'.join(self.Sequence[i:i+60] for i in range(0, len(self.Sequence), 60)) + "\n")


