# CAARS: Comparative Assembly and Annotation of RNA-Seq data

A method to introduce RNA-Seq data in existing multi-species multiple sequence alignments and reconstruct reliable phylogenies.

Any question or suggestion on the program can be addressed to: carine.rey@ens-lyon.org

For more information see the wiki page [https://github.com/CarineRey/caars/wiki](https://github.com/CarineRey/caars/wiki)

# Installation

## Test version

CAARS is available in a docker container available in DockerHub.

(If you don't have docker, you can find [here](https://docs.docker.com/linux/step_one/) installation instruction.)

You can get and run the container with this command:

```sh
cd /shared/directory
export SHARED_DIR=$PWD

# start the docker
docker run -t -i -e LOCAL_USER_ID=`id -u $USER` -e W_DIR=$SHARED_DIR -v $SHARED_DIR:$SHARED_DIR carinerey/caars bash
```

And you're done!


You **can** use:
 *  ``` -v $SHARED_DIR:$SHARED_DIR  ``` to share your working directory between your computer and the virtual environment in the docker container.
Indeed, CAARS builds links with absolute path which will be broken if you don't use the same directory tree.

 * ``` -e LOCAL_USER_ID=`id -u $USER` ``` to allow giving user rights on files created in the docker container.
 * ``` -e W_DIR=$SHARED_DIR ``` to set the working directory as the shared directory in the docker container.

CAARS can be called directly in the docker container terminal.

```sh
user_caars@1dbd8f2594cc:/shared/directory$ caars -help
```

All data used in the container must be placed in the $SHARED_DIRECTORY, indeed the container "sees" only directory tree from the shared directory. 

See the [Tutorial](https://github.com/CarineRey/caars/wiki/Tutorial) for more information.

If you have any problem do not hesitate to contact me (carine.rey@ens-lyon.org).


## Complete installation

See [INSTALL.md](INSTALL.md) or [INSTALL-HPC.md](INSTALL-HPC.md).

# Usage

Example:

```
caars  --outdir OUTPUT_DIR --sample-sheet sample_sheet.tsv --species-tree /home/user/data/species_tree.nw --alignment-dir GENE_FAMILIES_MSA_DIR --seq2sp-dir GENE_FAMILIES_SEQ2SP_DIR --np 2 --memory 5
```

CAARS will work in a directory named ```_caars``` build in the current directory.
This directory will contain all output files but under coded names.
At the end, CAARS will build symbolic links between the OUTPUT_DIR and the ```_caars`` directory.

To be sure not to loose your output file, copy the OUTPUT_DIR to a new directory by taking into account links.
For instance:
```
cp -rL OUTPUT_DIR OK_OUTPUT_DIR
```

See the [Tutorial](https://github.com/CarineRey/caars/wiki/Tutorial) for more usage.


```
caars

  caars

=== flags ===

  --alignment-dir PATH     Directory containing all gene family alignments
                           (Family_name.fa) in fasta format.
  --outdir PATH            Destination directory.
  --sample-sheet PATH      Sample sheet file.
  --seq2sp-dir PATH        Directory containing all link files
                           (Family_name.tsv). A line for each sequence and its
                           species spaced with tabulations.
  --species-tree ABSOLUTE  PATH Species tree file in nw format containing all
                           species. Warning absolute path is required.
  [--dag-graph PATH]       Write dag graph in an dot file (Can take a lot of
                           time)
  [--debug]                Get intermediary files (Default:false)
  [--html-report PATH]     Logs build events in an HTML report
  [--just-parse-input]     Parse input and exit. Recommended to check all input
                           files. (Default:false)
  [--memory INT]           Number of GB of system memory to use.(Default:1)
  [--mpast FLOAT]          Minimal percentage of alignment of a CAARS
                           sequences on its (non CAARS) closest sequence to be
                           kept in the final output
  [--no-reconcile]         Not run final Reconciliation step
  [--np INT]               Number of CPUs (at least 2). (Default:2)
  [--quiet]                Do not report progress. Default: off
  [--refinetree]           Refine topology during final reconciliation step
                           (Default:false)
  [-build-info]            print info about this build and exit
  [-version]               print the version of this build and exit
  [-help]                  print this help text and exit
                           (alias: -?)


```

### alignment-dir

A directory containing  fasta formated multiple sequence alignments (MSA).
Each filename must be composed of the name of the gene family and the extension ".fa"

### seq2sp-dir

A directory containing tabular link files between sequence and species.
Each filename must be composed of the name of the gene family and the extension ".tsv"

### sample-sheet

CAARS needs in input a sample sheet file.

This file is composed of 10 tabulated delimited columns with headers:
  * Sample ID: an unique identifier (3 capital letters is recommended)
  * Sample species name: the species of the sample
  * Reference species name: A reference species to annotate the sample.
  * Path for a single-end RNA-seq run : Single *fastq* or *fasta* file path (authorized extension: .fa, .fasta, .fq, .fastq)
  * Path for a double-end RNA-seq run : Left *fastq* or *fasta* file path
  * Path for a double-end RNA-seq run : Right *fastq* or *fasta* file path
  * Strand and type of the RNA-seq run : F,R,RF,FR,US,UP
  * Run the standard assembly on data : Yes or No
  * Path to a given assembly : *fasta* file path (must contain only cds)
  * Run apytram on data : Yes or No

Column order must be conserved and also the header line.

An example:

id	|species	|ref_species	|path_fastx_single	|path_fastx_left	|path_fastx_right	|orientation	|run_standard	|path_assembly	|run_apytram
---|---|---|---|---|---|---|---|---|---
CMA	|Mesocricetus_auratus	|Homo_sapiens	|-	|fastq/Mesocricetus_auratus.1.fq	|fastq/Mesocricetus_auratus.2.fq	|UP	|yes	|Trinity_assembly.AMH.fa	|yes
CMM	|Mus_musculus	|Homo_sapiens,Cavia_porcellus	|fasta/Mus_musculus.fa	| -	|-	|F	|yes	|-	|yes


### species-tree

**Absolute path** to the species tree (a topology can be enough) including all species

###  outdir
A directory path which will contain outputs


## Run CAARS on test datasets

Go to the [Tutorial](https://github.com/CarineRey/caars/wiki/Tutorial) page to have usage examples.

Two test datasets are also available in the source directory of CAARS.

* A very small dataset (1 min):
```sh
# In the CAARS directory (/opt/CAARS in the docker container) run:
make test
# or in the example directory (/opt/CAARS/example in the docker container)
bash Launch_CAARS.sh

# To remove outputs
make clean_test
```

* A small dataset (10 min):
```sh
# In the CAARS directory (/opt/CAARS in the docker container) run:
make test2
# or in the example directory (/opt/CAARS/example in the docker container):
bash Launch_CAARS2.sh

# To remove outputs
make clean_test2
```
