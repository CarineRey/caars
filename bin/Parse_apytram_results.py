#!/usr/bin/python
# coding: utf-8

import sys
import os
import logging
import re
import string

### Set up the logger
# create logger with 'spam_application'
logger = logging.getLogger('parse_apytram_input')
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


if len(sys.argv) != 3:
    logger.error("2 arguments are required")
    sys.exit(1)

ConfigFileName = sys.argv[1]
OutDirName = sys.argv[2]

if not os.path.isfile(ConfigFileName):
    logger.error("The first argument must be a file")
    sys.exit(1)


### Set up the output directory
if OutDirName:
    if os.path.isdir(OutDirName):
        logger.info("The output directory %s exists" %(OutDirName) )
    elif os.path.isfile(OutDirName):
        logger.error("The second argument must be a directory")
        sys.exit(1)
    elif OutDirName: # if OutDirName is not a empty string we create the directory
        logger.info("The output directory %s does not exist, it will be created" %(OutDirName))
        os.makedirs(OutDirName)
else:
    logger.error("The second argument must be a directory")
    sys.exit(1)


# A function to read and write a new fasta and a sp2seq file
def read_fasta_from_apytram(FastaPath2Sp_dic,
                            OutFastaFileName,
                            Sp2SeqFileName,
                            SeqId_dic
                            ):
    """Read a list of fasta and write a new fasta file with unique sequence names
    Build also a Sp2Seq lin file"""
    String_list = []
    SeqName_list = []
    
    for (InFastaFileName, Species) in FastaPath2Sp_dic.items():
		InFile = open(InFastaFileName,"r") 
		for line in InFile:
			if re.match(">",line):
				# This is a new sequence
				SeqName = "%s%s%s" %(SeqId_dic[Species]["SeqPrefix"],
				                     string.zfill(SeqId_dic[Species]["SeqNb"],
				                                  SeqId_dic[Species]["NbFigures"]
				                                  ),
				                     "_%s" %(Family))
				SeqId_dic[Species]["SeqNb"] += 1
				SeqName_list.append(SeqName)
				String_list.append( ">%s\n" %(SeqName))
			else:
				String_list.append(line)
		InFile.close()
    
    # Write all sequences in the output fasta file
    OutFile = open(OutFastaFileName,"w")
    OutFile.write("".join(String_list))
    OutFile.close()
    
    #Build and write the Sp2SeqFile
    Sp2SeqFile = open(Sp2SeqFileName,"w")
    Sp2Seq_String = "%s:" %(Species) + ("\n%s:" %(Species)).join(SeqName_list)
    Sp2SeqFile.write(Sp2Seq_String)
    Sp2SeqFile.close()
    return(SeqId_dic)


### Read the config file:
## Sp\tFam\tdir_path
ConfigFile = open(ConfigFileName,"r")
SeqId_dic = {}
FastaPath2SpPerFam_dic = {}
ApytramPrefix = "apytram"
NbFigures = 10
for line in ConfigFile:
    (Species, Family, ApytramDir) = line.strip().split("\t")
    InFastaFileName = "%s/%s.fasta" %(ApytramDir,ApytramPrefix)
    if os.path.isfile(InFastaFileName):
		FastaPath2SpPerFam_dic.setdefault(Family,{})
		FastaPath2SpPerFam_dic[Family][InFastaFileName] = Species
		SeqId_dic.setdefault(Species, {"SeqPrefix" : "AP%s0" %(Species[0:3].upper()),
									   "SeqNb" : 1, "NbFigures" : NbFigures})

for Family in FastaPath2SpPerFam_dic.keys():
	OutFastaFileName = "%s/apytram.%s.fa" %(OutDirName,Family)
	Sp2SeqFileName = "%s/apytram.%s.sp2seq.txt" %(OutDirName,Family)
	SeqId_dic = read_fasta_from_apytram(FastaPath2SpPerFam_dic[Family],
                            OutFastaFileName,
                            Sp2SeqFileName,
                            SeqId_dic)
                            
sys.exit(0)

    
    
