#!/bin/bash

CURRENT_DIR=$PWD
WORKING_DIR=$CURRENT_DIR/working_dir
OUTPUT_DIR=$CURRENT_DIR/output_dir
DATA_DIR=$CURRENT_DIR/data


echo "START:"
date
START=$(date +%s)

mkdir -p $WORKING_DIR
cd $WORKING_DIR

amalgam_app.byte  --outdir $OUTPUT_DIR --sample-sheet $DATA_DIR/sample_sheet.tsv --species-tree $DATA_DIR/species_tree.nw --alignment-dir $DATA_DIR/gene_fams/ --seq2sp-dir $DATA_DIR/sp2seq_links/ --np 2 --memory 5 --mpast 50 

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "It took $DIFF seconds"
echo "END:"
date

