
import os
import sys
import logging
import subprocess


class Makeblastdb:
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
        ExitCode = 1
        command = ["makeblastdb","-in",self.InputFile,"-out", os.path.abspath(self.OutputFiles),
                   "-dbtype", self.Dbtype ]
        if self.IndexedDatabase:
            command.append("-parse_seqids")
        
        self.logger.debug(" ".join(command))
        p = subprocess.Popen(command,
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE)
        out, err = p.communicate()
        if err:
            self.logger.error(err)
        return (out,err)
    

class Blast:
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
   
    def launch(self,OutputFile):
        ExitCode = 1
        if self.Program in ["blastn","blastx","tblastn","tblastx"]:
            command = [self.Program,"-db", self.Database, "-query" , self.QueryFile,
                      "-evalue", str(self.Evalue), "-outfmt", str(self.OutFormat),
                      "-out", OutputFile,
                      "-max_target_seqs", str(self.max_target_seqs),
                      "-num_threads", str(self.Threads)]
            
            if self.perc_identity:
                command.extend(["-perc_identity" , str(self.perc_identity)])
            if self.Task:
                command.extend(["-task",self.Task])
                
            self.logger.debug(" ".join(command))
            p = subprocess.Popen(command,
                                           stdout=subprocess.PIPE,
                                           stderr=subprocess.PIPE)
            out, err = p.communicate()
            if err:
                self.logger.error(err)
            return (out,err)

        else:
                self.logger.error("%s not in [blastn,blastx,tblastn,tblastx]" % self.Program)
                return ("","%s not in [blastn,blastx,tblastn,tblastx]" % self.Program)
    
class Blastdbcmd:
    """Define a object to launch blastdbcmd on a local database"""
    def __init__(self, Database, SequenceNamesFile, OutputFile):
        self.logger = logging.getLogger('main.lib.BlastPlus.Blastdbcmd')
        self.InputFile = SequenceNamesFile
        self.Database = Database
        self.OutputFile = OutputFile
        self.Dbtype = "nucl"

    def launch(self):
        command = ["blastdbcmd","-db",self.Database,"-entry_batch", self.InputFile,
                   "-dbtype", self.Dbtype]
                   
        if self.OutputFile:
			command.extend(["-out" , self.OutputFile ])
        self.logger.debug(" ".join(command))
        p = subprocess.Popen(command,
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE)
        out, err = p.communicate()
        if err:
            self.logger.error(err)
        return (out,err)
    
    def is_database(self):
        Out = False
        command = ["blastdbcmd","-db",self.Database,"-info"]
        self.logger.debug(" ".join(command))
        p = subprocess.Popen(command,
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE)
        out, err = p.communicate()
        if err:
            self.logger.error(err)
        else:
       # If no error the database exist
            Out = True             
        return Out
