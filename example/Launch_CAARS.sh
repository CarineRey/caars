#!/bin/bash

set -euo pipefail +o nounset

CURRENT_DIR=$PWD
WORKING_DIR=$CURRENT_DIR/working_dir
OUTPUT_DIR=$CURRENT_DIR/output_dir
DATA_DIR=$CURRENT_DIR/data


echo "START:"
date
START=$(date +%s)

mkdir -p $WORKING_DIR
cd $WORKING_DIR

if (( $# == 0 ))
then
    echo "no arguments"
    use_docker=""
else
    if [[ $1 == "docker" ]]
        then
        use_docker="--use-docker"
        else
        use_docker=""
        fi
fi

../../caars  --outdir $OUTPUT_DIR/justparseinput --sample-sheet $DATA_DIR/sample_sheet.tsv --species-tree $DATA_DIR/species_tree.nw --alignment-dir $DATA_DIR/gene_fams/ --seq2sp-dir $DATA_DIR/sp2seq_links/ --np 2 --memory 5 --mpast 50  $use_docker --family-subset $DATA_DIR/subset_fam.txt --just-parse-input
../../caars  --outdir $OUTPUT_DIR/all --sample-sheet $DATA_DIR/sample_sheet.tsv --species-tree $DATA_DIR/species_tree.nw --alignment-dir $DATA_DIR/gene_fams/ --seq2sp-dir $DATA_DIR/sp2seq_links/ --np 2 --memory 5 --mpast 50  $use_docker --family-subset $DATA_DIR/subset_fam.txt

files="$OUTPUT_DIR/justparseinput/UsableFamilies.txt
$OUTPUT_DIR/all/assembly_results_only_seq/CAARS_sequences.seq2sp2fam.txt"

for file in $files
do
if [ -s "$file" ]
then
    l="`wc -l < $file`"
    if [ "$l" -gt "1" ]
    then
        echo "$file has some data."
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

