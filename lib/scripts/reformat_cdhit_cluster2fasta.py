#!/usr/bin/python
# coding: utf-8

# File: parse_clustr.py
# Created by: Carine Rey
# Created on: March 2017


import itertools,time,logging,sys

### Set up the logger
# create logger with 'spam_application'
logger = logging.getLogger('parse_clustr')
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
logger.debug(" ".join(sys.argv))

if len(sys.argv) != 3:
    logger.error("2 arguments are required (clustr file, output file)")
    sys.exit(1)

clstr_filename = sys.argv[1]
output_filename = sys.argv[2]

logger.debug("Clustr file: %s", clstr_filename)
logger.debug("Output file: %s", output_filename)

start_time = time.time()

with open(output_filename, "w") as output_file:
    with open(clstr_filename, "r") as cluster_file:
        cluster_groups = (x[1] for x in itertools.groupby(cluster_file, key=lambda line: line[0] == '>'))
        print("%s seconds" %( str(time.time() - start_time)))
        print "write ouptut"
        for cluster in cluster_groups:
            seqs = [seq.split('>')[1].split('...')[0] for seq in cluster_groups.next()]
            if len(seqs) > 1:
                output_file.write(">%s\n%s\n" %(seqs[0], ";".join(seqs[1:])))

    print("--- %s seconds ---" %( str(time.time() - start_time)))
