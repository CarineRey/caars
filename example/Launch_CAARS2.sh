#!/bin/bash

set -euo pipefail

CURRENT_DIR=$PWD
WORKING_DIR=$CURRENT_DIR/working2_dir
OUTPUT_DIR=$CURRENT_DIR/output2_dir
DATA_DIR=$CURRENT_DIR/data

if [[ $1 == "docker" ]]
then
use_docker="--use-docker"
else
use_docker=""
fi


echo "START:"
date
START=$(date +%s)

mkdir -p $WORKING_DIR
cd $WORKING_DIR

echo "_no_filter"
../../caars  --outdir $OUTPUT_DIR"_no_filter"         $use_docker --sample-sheet $DATA_DIR/sample_sheet.tsv --species-tree $DATA_DIR/species_tree.nw --alignment-dir $DATA_DIR/gene_fams/ --seq2sp-dir $DATA_DIR/sp2seq_links/ --np 3 --memory 5

echo "_no_reconcilation"
../../caars  --outdir $OUTPUT_DIR"_no_reconcilation"  $use_docker --sample-sheet $DATA_DIR/sample_sheet.tsv --species-tree $DATA_DIR/species_tree.nw --alignment-dir $DATA_DIR/gene_fams/ --seq2sp-dir $DATA_DIR/sp2seq_links/ --np 3 --memory 5 --no-reconcile

echo "_refinetree"
../../caars  --outdir $OUTPUT_DIR"_refinetree"        $use_docker --sample-sheet $DATA_DIR/sample_sheet.tsv --species-tree $DATA_DIR/species_tree.nw --alignment-dir $DATA_DIR/gene_fams/ --seq2sp-dir $DATA_DIR/sp2seq_links/ --np 3 --memory 5 --refinetree

echo "_merge"
../../caars  --outdir $OUTPUT_DIR"_merge"             $use_docker --sample-sheet $DATA_DIR/sample_sheet.tsv --species-tree $DATA_DIR/species_tree.nw --alignment-dir $DATA_DIR/gene_fams/ --seq2sp-dir $DATA_DIR/sp2seq_links/ --np 3 --memory 5 --mpast 50

echo "_length"
../../caars  --outdir $OUTPUT_DIR"_length"            $use_docker --sample-sheet $DATA_DIR/sample_sheet.tsv --species-tree $DATA_DIR/species_tree.nw --alignment-dir $DATA_DIR/gene_fams/ --seq2sp-dir $DATA_DIR/sp2seq_links/ --np 3 --memory 5 --mpast 50 --merge-criterion length --no-reconcile

echo "_length_complete"
../../caars  --outdir $OUTPUT_DIR"_length_complete"   $use_docker --sample-sheet $DATA_DIR/sample_sheet.tsv --species-tree $DATA_DIR/species_tree.nw --alignment-dir $DATA_DIR/gene_fams/ --seq2sp-dir $DATA_DIR/sp2seq_links/ --np 3 --memory 5 --mpast 50 --merge-criterion length_complete --no-reconcile

echo "_big"
../../caars  --outdir $OUTPUT_DIR"_big"               $use_docker --sample-sheet $DATA_DIR/sample_sheet2.tsv --species-tree $DATA_DIR/species_tree.nw --alignment-dir $DATA_DIR/gene_fams/ --seq2sp-dir $DATA_DIR/sp2seq_links/ --np 3 --memory 5 --mpast 50

SUFFIX_l="_no_filter
_no_reconcilation
_merge
_refinetree
_length
_big
_length_complete"


for SUFFIX in $SUFFIX_l
do
file=$OUTPUT_DIR$SUFFIX/assembly_results_only_seq/CAARS_sequences.seq2sp2fam.txt

if [ -s "$file" ]
then
    l="`wc -l < $file`"
    if [ "$l" -gt "1" ]
    then
        echo "$file has some data."
        echo "... $SUFFIX => ... OK !"
    else
        echo "$file has no data."
        exit 1
    fi
else
    echo "$file is empty."
    exit 1
fi
done

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "It took $DIFF seconds"
echo "END:"
date

