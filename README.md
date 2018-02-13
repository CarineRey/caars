# CAARS: Comparative Assembly and Annotation of RNA-Seq data

[![Build Status](https://travis-ci.org/CarineRey/caars.svg?branch=master)](https://travis-ci.org/CarineRey/caars)

A method to introduce RNA-Seq data in existing multi-species multiple sequence alignments and reconstruct reliable phylogenies.

Any question or suggestion on the program can be addressed to: carine.rey [at] ens-lyon.org

For more information see the wiki page [https://github.com/CarineRey/caars/wiki](https://github.com/CarineRey/caars/wiki)

or the [Tutorial](https://github.com/CarineRey/caars/wiki/Tutorial).

# Installation

## No installation required with Docker and docker container usage

CAARS is available in a docker image available on DockerHub.

(If you don't have docker, you can find [here](https://docs.docker.com/install/) installation instruction.)

You can get the docker image and run CAARS through a docker container with these commands:

```{sh}
## 1/ Get (or update) the caars docker image

docker pull carinerey/caars

## 2/ On your machine, go to the working directory

cd /shared/directory
export SHARED_DIR=$PWD


## 3/ Use caars through the docker container

# Define a function to simplify the docker container usage

function docker_caars_cmd { echo "docker run --rm -e LOCAL_USER_ID=`id -u $USER` -v $SHARED_DIR:$SHARED_DIR -w `pwd` carinerey/caars "
}

# "-e LOCAL_USER_ID=`id -u $USER`"  will allow giving user rights on files created in the docker container.
# "-v $SHARED_DIR=$SHARED_DIR"      will allow sharing the $SHARED_DIR directory between your computer and the virtual environment in the docker container.
# "-w $SHARED_DIR"                  will set the working directory as the $SHARED_DIR directory in the docker container.

# Use the function to run caars

`docker_caars_cmd` caars [options]

# For example to get the help message:

`docker_caars_cmd` caars  -h

```

And you're done!


All data used in the container must be placed in the $SHARED_DIRECTORY, indeed the container detects/recognizes only directory tree from the shared directory. 

If you have any problem do not hesitate to contact me (carine.rey [at] ens-lyon.org).


## Local installation

See [INSTALL.md](INSTALL.md) or [INSTALL-HPC.md](INSTALL-HPC.md).


# Usage


## Basic command line:

```
caars  --outdir OUTPUT_DIR --sample-sheet sample_sheet.tsv --species-tree /home/user/data/species_tree.nw --alignment-dir GENE_FAMILIES_MSA_DIR --seq2sp-dir GENE_FAMILIES_SEQ2SP_DIR --np 2 --memory 5
```

CAARS will work in a directory named ```_caars``` build in the current directory.
This directory will contain all output files but under coded names.
At the end, CAARS will build symbolic links between the OUTPUT_DIR and the ```_caars``` directory.

To be sure not to loose your output files, copy the OUTPUT_DIR to a new directory by taking into account links.
For instance:
```
cp -rL OUTPUT_DIR OK_OUTPUT_DIR
```

## All options:

```
caars

  caars 

=== flags ===

  --alignment-dir PATH     Directory containing all gene family alignments
                           (Family_name.fa) in fasta format.
  --outdir PATH            Destination directory.
  --sample-sheet PATH      sample sheet file.
  --seq2sp-dir PATH        Directory containing all linked files
                           (Family_name.tsv). Each line contains a sequence name
                           and its species spaced by a tabulation.
  --species-tree ABSOLUTE  PATH Species tree file in nw format containing all
                           species. Warning absolute path is required.
  [--dag-graph PATH]       Write dag graph in a dot file (Can take a lot of
                           time)
  [--debug]                Get intermediary files (Default:false)
  [--family-subset PATH]   A file containing a subset of families to use.
                           Default: off
  [--html-report PATH]     Logs build events in an HTML report
  [--just-parse-input]     Parse input and exit. Recommended to check all input
                           files. (Default:false)
  [--memory INT]           Number of GB of system memory to use.(Default:1)
  [--merge-criterion STR]  Merge criterion during redundancy removing. It must
                           be “length“ or “length_complete” or
                           “merge”. “length” means the longest sequence
                           is selected. “length.complete” : means the
                           largest number of complete sites (no gaps).
                           “merge” means that the set of monophyletic
                           sequences is used to build one long “chimera”
                           sequence corresponding to the merging of them.
  [--mpast FLOAT]          Minimal percentage of alignment of a caars sequence
                           on its (non Caars) closest sequence to be kept in the
                           final output
  [--no-reconcile]         Not run final Reconciliation step
  [--np INT]               Number of CPUs (at least 2). (Default:2)
  [--quiet]                Do not report progress. Default: off
  [--refinetree]           Refine topology during final Reconciliation step
                           (Default:false)
  [--use-docker]           Use docker in caars. Default: off
  [-build-info]            print info about this build and exit
  [-version]               print the version of this build and exit
  [-help]                  print this help text and exit
```


## Required options

### alignment-dir

A directory containing fasta formated multiple sequence alignments (MSA).
Each filename must be composed of the name of the gene family and the extension ".fa".
**Warning:** Alignments must only contain "ATGCNUWSMKRYBDHV-" characters.

### seq2sp-dir

A directory containing tabular files including links between sequences and species.
Each filename must be composed of the name of the gene family and the extension ".tsv".

### sample-sheet

CAARS needs a sample sheet file.

This file is composed of 10 tabulation-delimited columns with headers:
  * Sample ID: a unique identifier (3 capital letters is recommended)
  * Sample species name: the species of the sample
  * Group ID: the group of the sample
  * Reference species name: Reference species used to annotate the sample.
  * Path for a single-end RNA-seq run : Single *fastq* or *fasta* file path (authorized extension: .fa, .fasta, .fq, .fastq)
  * Path for a paired-end RNA-seq run : Left *fastq* or *fasta* file path
  * Path for a paired-end RNA-seq run : Right *fastq* or *fasta* file path
  * Strand and type of the RNA-seq run : F,R,RF,FR,US or UP
  * Run the standard assembly on data : Yes or No
  * Path to a given assembly : *fasta* file path (must contain only cds)
  * Run the assisted assembly (apytram) on data : Yes or No

Column order must be conserved and also the header line.

An example:

id	|species	|group_id	|ref_species	|path_fastx_single	|path_fastx_left	|path_fastx_right	|orientation	|run_standard	|path_assembly	|run_apytram
---|---|---|---|---|---|---|---|---|---|---
CMA	|Mesocricetus_auratus	|group1	|Homo_sapiens	|-	|fastq/Mesocricetus_auratus.1.fq	|fastq/Mesocricetus_auratus.2.fq	|UP	|yes	|Trinity_assembly.AMH.fa	|yes
CMM	|Mus_musculus	|group2	|Homo_sapiens,Cavia_porcellus	|fasta/Mus_musculus.fa	| -	|-	|F	|yes	|-	|yes


### species-tree

**Absolute path** to the species tree (a topology can be enough) including all species

### outdir

A directory path which will contain outputs



# Test dataset and Tutorial

## Run CAARS on a test dataset

Go to the [Tutorial](https://github.com/CarineRey/caars/wiki/Tutorial) page to have usage examples.

A test dataset is also available in the source directory of CAARS.

* A small dataset (10 min):
```{sh}
# In the CAARS directory (/opt/CAARS in the docker container) run:
make test
# or in the example directory (/opt/CAARS/example in the docker container)
bash Launch_CAARS.sh

# To remove outputs
make clean_test
```

