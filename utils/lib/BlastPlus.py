# File: BlastPlus.py
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
import logging
import subprocess


class Makeblastdb(object):
    """Define an object to create a local database"""
    def __init__(self, InputFile, OutputFiles):
        self.logger = logging.getLogger('main.lib.BlastPlus.Makeblastdb')
        self.Title = ""
        self.InputFile = InputFile
        self.OutputFiles = OutputFiles
        self.LogFile = ""
        self.Dbtype = "nucl"
        self.IndexedDatabase = True

    def launch(self):
        command = ["makeblastdb", "-in", self.InputFile,
                   "-out", os.path.abspath(self.OutputFiles),
                   "-dbtype", self.Dbtype]
        if self.IndexedDatabase:
            command.append("-parse_seqids")

        self.logger.debug(" ".join(command))
        p = subprocess.Popen(command,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        (out, err) = p.communicate()
        if err:
            self.logger.error(err)
        return (out, err)


class Blast(object):
    """Define a object to lauch a blast on a local database"""
    def __init__(self, Program, Database, QueryFile):
        self.logger = logging.getLogger('main.lib.BlastPlus.Blast')
        self.Program = Program
        self.Database = Database
        self.QueryFile = QueryFile
        self.OutputFile = ""
        self.Threads = 1
        self.OutFormat = 8
        self.Evalue = 10
        self.max_target_seqs = 1000000000
        self.perc_identity = 0
        self.Task = ""

    def launch(self, OutputFile):
        if self.Program in ["blastn", "blastx", "tblastn", "tblastx"]:
            command = [self.Program, "-db", self.Database,
                      "-query", self.QueryFile,
                      "-evalue", str(self.Evalue),
                      "-outfmt", str(self.OutFormat),
                      "-out", OutputFile,
                      "-max_target_seqs", str(self.max_target_seqs),
                      "-num_threads", str(self.Threads)]

            if self.perc_identity:
                command.extend(
                ["-perc_identity", str(self.perc_identity)]
                )
            if self.Task:
                command.extend(["-task", self.Task])

            self.logger.debug(" ".join(command))
            p = subprocess.Popen(command,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
            (out, err) = p.communicate()
            if err:
                self.logger.error(err)
            return (out, err)

        else:
            self.logger.error(
            "%s not in [blastn,blastx,tblastn,tblastx]",
            self.Program
            )
            return ("", "%s not in [blastn,blastx,tblastn,tblastx]" %self.Program)

class Blastdbcmd(object):
    """Define a object to launch blastdbcmd on a local database"""
    def __init__(self, Database, SequenceNamesFile, OutputFile):
        self.logger = logging.getLogger('main.lib.BlastPlus.Blastdbcmd')
        self.InputFile = SequenceNamesFile
        self.Database = Database
        self.OutputFile = OutputFile
        self.Dbtype = "nucl"

    def launch(self):
        command = ["blastdbcmd", "-db", self.Database,
                   "-entry_batch", self.InputFile,
                   "-dbtype", self.Dbtype]

        if self.OutputFile:
            command.extend(["-out", self.OutputFile])
        self.logger.debug(" ".join(command))
        p = subprocess.Popen(command,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        (out, err) = p.communicate()
        if err:
            self.logger.error(err)
        return (out, err)

    def is_database(self):
        Out = False
        command = ["blastdbcmd", "-db", self.Database, "-info"]
        self.logger.debug(" ".join(command))
        p = subprocess.Popen(command,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        (_, err) = p.communicate()
        if err:
            self.logger.info(err)
        else:
       # If no error the database exist
            Out = True
        return Out
