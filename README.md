# Amalgam

A method to introduce RNA-seq data in existing multi-species multi asequences alignments and reconstruct reliable phylogenies

[ In development ]

Any question or suggestion on the program can be addressed to: carine.rey@ens-lyon.fr

For more information gto to the wiki page [https://github.com/CarineRey/amalgam/wiki](https://github.com/CarineRey/amalgam/wiki)

# Installation

## Test version

amalgam is available in a docker container available in DockerHub.

(If you don't have docker, you can find [here](https://docs.docker.com/linux/step_one/) installation instruction.)

You can get and run the container with this command:

```sh
docker run -t -i -v $pwd:$pwd carinerey/amalgam
```

And you're done!

You **must** use the -v option to share your working directory between your computer and the virtual environment in the docker container. Indeed, amalgam builds links with absolute path which will be break if you don't use the same directory tree.

amalgam can be call directly in the docker container terminal.

Warning, the container is in development, if you have any problem don't hesitate to contact me (carine.rey@ens-lyon.fr).

## complete installation

See [INSTALL.md](INSTALL.md) or [INSTALL-HPC.md](INSTALL-HPC.md).

# Usage

Example:

```
amalgam_app.byte  --outdir OUTPUT_DIR --sample-sheet sample_sheet.tsv --species-tree /home/user/data/species_tree.nw --alignment-dir GENE_FAMILIES_MSA_DIR --seq2sp-dir GENE_FAMILIES_SEQ2SP_DIR --np 2 --memory 5
```

Amalgam will work in a directory named ```_bistro``` build in the current directory.
This directory will contain all output files but under coded names.
At the end, amalgam will build symbolic links between the OUTPUT_DIR and the ```_bistro``` directory.

To be sure to not lost your output file, copy the OUTPUT_DIR to a new directory by taking into account links.
For instance:
```
cp -rL OUTPUT_DIR OK_OUTPUT_DIR
```

```
amalgam_app.byte

=== flags ===

  --alignment-dir PATH        Directory containing all gene family alignments
                              (Family_name.fa) in fasta format.
  --outdir PATH               Destination directory.
  --sample-sheet PATH         sample sheet file.
  --seq2sp-dir PATH           Directory containing all link files
                              (Family_name.tsv). A line for each sequence and
                              its species spaced by a tabulation.
  --species-tree ABSOLUTE     PATH Species tree file in nw format containing all
                              species. Warning absolute path is required.
  [--dag-graph PATH]          Write dag graph in an dot file (Can take a lot of
                              time)
  [--debug Get]               intermediary files (Default:false)
  [--html-report PATH]        Logs build events in an HTML report
  [--just-parse-input Parse]  input and exit. Recommended to check all input
                              files. (Default:false)
  [--memory INT]              Number of GB of system memory to use.(Default:1)
  [--mpast FLOAT]             Minimal percentage of alignment of an Amalgam
                              sequences on its (non amalgam) closest sequence to
                              be kept in the final output
  [--no-reconcile Not]        run final Reconciliation step
  [--np INT]                  Number of CPUs (at least 2). (Default:2)
  [--quiet Do]                not report progress. Default: off
  [--refinetree Refine]       topology during final Reconciliation step
                              (Default:false)
  [-build-info]               print info about this build and exit
  [-version]                  print the version of this build and exit
  [-help]                     print this help text and exit
                              (alias: -?)

```

### alignment-dir

A directory containing  fasta formated multi sequence alignments (MSA).
Each filename must be composed of the name of the gene family and the extension ".fa"

### seq2sp-dir

A directory containing tabular link files between sequence and species.
Each filename must be composed of the name of the gene family and the extension ".tsv"

### sample-sheet

amalgam needs in input a sample sheet file.

This file is composed of 10 tabulated delimited columns with headers:
  * Sample ID: an unique identifiant (3 capital letters is recommended)
  * Sample species name: the species of the sample
  * Reference species name: A reference species to annotate the sample.
  * Path for a single-end RNA-seq run : Single *fastq* file path
  * Path for a double-end RNA-seq run : Left *fastq* file path
  * Path for a double-end RNA-seq run : Right *fastq* file path
  * Strand and type of the RNA-seq run : F,R,RF,FR,US,UP
  * Run the standard on data : Yes or No
  * Path to a given trinity assembly : fasta file path
  * Run apytram on data : Yes or No

Column order must be conserved and also the header line.

An example:

id	|species	|ref_species	|path_fastq_single	|path_fastq_left	|path_fastq_right	|orientation	|run_trinity	|path_assembly	|run_apytram
---|---|---|---|---|---|---|---|---|---
AMH	|Mesocricetus_auratus	|Homo_sapiens	|-	|fastq/Mesocricetus_auratus.1.fq	|fastq/Mesocricetus_auratus.2.fq	|UP	|yes	|Trinity_assembly.AMH.fa	|yes
AMM	|Mus_musculus	|Homo_sapiens,Cavia_porcellus	|fastq/Mus_musculus.fq	| -	|-	|F	|yes	|-	|yes


### species-tree

**Absolute path** to the species tree (a topology can be enough) including all species

###  outdir
A directory path which will contain outputs


## Run amalgam on test datasets

Go to the [Tutorial](https://github.com/CarineRey/amalgam/wiki/Tutorial) page to have usage examples.

Two test datasets are also available in the source directory of amalgam.

* A very little dataset (1 min):
```
# In the amalgam directory run:
make test
# or in the example directory
bash Launch_amalgam.sh

# To remove outputs
make clean_test
```

* A  little dataset (10 min):
```
# In the amalgam directory run:
make test2
# or in the example directory
bash Launch_amalgam2.sh

# To remove outputs
make clean_test2
```
