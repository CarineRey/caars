# File: Aligner.py
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


class Exonerate(object):
    """Define an object to launch Exonerate"""
    def __init__(self, TargetFile, QueryFile):
        self.TargetFile = TargetFile
        self.QueryFile = QueryFile
        self.Model = ""
        self.Ryo = ""
        self.ShowAlignment = "F"
        self.ShowVulgar = "F"
        self.Bestn = 0
        self.Identity = 0
        self.Exhaustive = False

    def get_output(self):
        #Out = "Command line: []\nHostname:\n-- completed exonerate analysis"
        command = ["exonerate", "-t", self.TargetFile,
                   "-q", self.QueryFile,
                   "--showalignment", self.ShowAlignment,
                   "--showvulgar", self.ShowVulgar]

        if self.Model:
            command.extend(["--model", str(self.Model)])
        if self.Bestn > 0:
            command.extend(["--bestn", str(self.Bestn)])
        if self.Identity > 0:
            command.extend(["--percent", str(self.Identity)])
        if self.Ryo:
            command.extend(["--ryo", str(self.Ryo)])
        if self.Exhaustive:
            command.extend(["--exhaustive", "T"])

        try:
            Out = subprocess.check_output(command,
                                         stderr=open("/dev/null", "w"))
        except:
            os.system("echo Unexpected error when we launch Exonerate:\n")
            print " ".join(command)

        # Remove Exonerate default lines
        #(First, second (which contains Hostname) and last):
        List = Out.strip().split("\n")
        Start = 0
        while not "Hostname" in List[Start]:
            Start += 1

        Out = "\n".join(List[Start+1:-1])

        return Out


class Mafft(object):
    """Define an object to launch Mafft"""
    def __init__(self, InputFile):
        self.logger = logging.getLogger("main.lib.mafft")
        self.logger.info('creating an instance of Mafft')
        self.InputFile = InputFile
        self.OutputFile = ""
        self.AddOption = False
        self.AdjustdirectionOption = False
        self.AutoOption = False
        self.QuietOption = False

    def launch(self, output=""):
        command = ["mafft"]

        if self.AdjustdirectionOption:
            command.append("--adjustdirection")
        if self.AutoOption:
            command.append("--auto")
        if self.AddOption:
            if os.path.isfile(self.AddOption):
                command.extend(["--add", self.AddOption])
        if self.QuietOption:
            command.append("--quiet")

        if output or self.OutputFile:
            if output:
                self.OutputFile = output
            command.extend(["--out", self.OutputFile])

        command.append(self.InputFile)
        self.logger.debug(" ".join(command))

        p = subprocess.Popen(command,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        (out, err) = p.communicate()
        if err:
            self.logger.error(err)

        return (out, err)
