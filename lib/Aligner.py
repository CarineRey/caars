import os
import sys
import subprocess

## TO DO: check if Exonerate is installed and is path is Exonerate


class Exonerate:
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
        Out = "Command line: []\nHostname:\n-- completed exonerate analysis"
        command = ["exonerate","-t", self.TargetFile,"-q", self.QueryFile,
                   "--showalignment", self.ShowAlignment,
                   "--showvulgar", self.ShowVulgar]
        
        if self.Model:
            command.extend(["--model",str(self.Model)])
        if self.Bestn > 0:
            command.extend(["--bestn",str(self.Bestn)])
        if self.Identity > 0:
            command.extend(["--percent",str(self.Identity)])
        if self.Ryo:
            command.extend(["--ryo",str(self.Ryo)])
        if self.Exhaustive:
            command.extend(["--exhaustive","T"])

        try:
            Out = subprocess.check_output(command,
                                         stderr=open("/dev/null", "w"))
        except:
            os.system("echo Unexpected error when we launch Exonerate:\n")
            print " ".join(command)

        # Remove Exonerate default lines (First, second (which contains Hostname) and last):
        List = Out.strip().split("\n")
        Start = 0
        while not "Hostname" in List[Start]:
            Start += 1
            
        Out = "\n".join(List[Start+1:-1])
        
        return Out


class Mafft:
    """Define an object to launch Mafft"""
    def __init__(self, InputFile):
        self.InputFile = InputFile
        self.AddOption = False
        self.AdjustdirectionOption = False
        self.AutoOption = False
        self.QuietOption = False
        
    def get_output(self):
        Out = ""
        command = ["mafft"]

        if self.AdjustdirectionOption:
            command.append("--adjustdirection")
        if self.AutoOption:
            command.append("--auto")
        if self.AddOption:
            if os.path.isfile(self.AddOption):
                command.extend(["--add",self.AddOption])
        if self.QuietOption:
            command.append("--quiet")
        
        command.append(self.InputFile)
        try:
            Out = subprocess.check_output(command,
                                          stderr=open("/dev/null", "w"))
        except:
            os.system("echo Unexpected error when we launch Mafft:\n")
            print " ".join(command)
            
        return Out


    