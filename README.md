# amalgam
A method to introduce RNA-seq data in existing multi-species multi asequences alignments and reconstruct reliable phylogenies

[ In development ]

Any question or suggestion on the program can be addressed to: carine.rey@ens-lyon.fr

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

  --alignment-dir PATH     Directory containing all gene family alignments
                           (Family_name.fa) in fasta format.
  --outdir PATH            Destination directory.
  --sample-sheet PATH      sample sheet file.
  --seq2sp-dir PATH        Directory containing all link files
                           (Family_name.tsv). A line for each sequence and its
                           species spaced by a tabulation.
  --species-tree ABSOLUTE  PATH Species tree file in nw format containing all
                           species. Warning absolute path is required.
  [--memory INT]           Number of GB of system memory to use.(Default:1)
  [--np INT]               Number of CPUs (at least 2). (Default:2)
  [-build-info]            print info about this build and exit
  [-version]               print the version of this build and exit
  [-help]                  print this help text and exit

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
  * Path for a single-end RNA-seq run : Single fastq file path
  * Path for a double-end RNA-seq run : Left fastq file path
  * Path for a double-end RNA-seq run : Right fastq file path
  * Strand and type of the RNA-seq run : F,R,RF,FR,US,UP
  * Run trinity on data : Yes or No
  * Path to a given trinity assembly : fasta file path
  * Run apytram on data : Yes or No

 Column order must be conserved and also the header line.

An example:

id	|species	|ref_species	|path_fastq_single	|path_fastq_left	|path_fastq_right	|orientation	|run_trinity	|path_assembly	|run_apytram
---|---|---|---|---|---|---|---|---|---
AMH	|Mesocricetus_auratus	|Homo_sapiens	|-	|fastq/Mesocricetus_auratus.1.fq	|fastq/Mesocricetus_auratus.2.fq	|UP	|yes	|Trinity_assembly.AMH.fa	|yes
AMM	|Mus_musculus	|Homo_sapiens	|fastq/Mus_musculus.fq	| -	|-	|F	|yes	|-	|yes


### species-tree

**Absolute path** to the species tree (a topology can be enough) including all species

###  outdir
A directory path which will contain outputs


## Run amalgam on test datasets

Two test datasets is available in the source directory of amalgam.

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



[To be removed]

### Dependencies

* Transdecoder >= 3.0.1

* SRAToolKit >= 2.8.1-2

* FastTree >= 2.1.7 (Warning: The executable must be fasttree and not FastTree)

* Python 2.7 (with pip and setuptools)
    * PyQt4
    * SciPy
    * MySQLdb
    * lxml
    * ete2
    * ete3
    * profileNJ
    * pandas

* Trinity >=2.3
    * Java >= 1.8 (OpenJRE works)
    * Bowtie >= 2 (tested with 2.2.9)

* phyldog >= 1.1.0
    * libPLL >= 1.0.2 sequential
    * boost 1.55 <= . >= 1.49
    * bpp >= 2.2.0 (Bio++)

* PhyloMerge (0.2 from 2017/01)
    * bpp >= 2.2.0 (Bio++)

* apytram >= 1.0
    * exonerate >= 2.2.0
    * mafft >=7.1
    * blast+ >= 2.6
    * python = 2.7
    * Trinity >=2.3

* OCaml >= 4.03.0
    * bistro
        * oasis
        * solvuu-build
        * ocamlgraph


```
# to get sources
git clone https://github.com/CarineRey/amalgam


# to install amalgam
cd amalgam
make

# to test all dependencies
make test
```
