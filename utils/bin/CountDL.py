#!/usr/bin/python
# coding: utf-8

import argparse
import time
import sys
import logging
import glob


from ete3 import Tree

##########
# inputs #
##########
start_time = time.time()

### Option defining
parser = argparse.ArgumentParser(prog="CountDL.py",
                                 description='')
parser.add_argument('--version', action='version', version='%(prog)s 1.0.0')

parser.add_argument('-sp_tree', type=str, help="species tree", required=True)
parser.add_argument('-o', type=str, help="species tree with the number of DL", required=True)
parser.add_argument('-rec_trees_dir', type=str, help="A file containing reconciled tree paths", required=True)
parser.add_argument('--debug', action="store_true", help="debug", default=False)

### Option parsing
args = parser.parse_args()

### Set up the logger
# create logger
logger = logging.getLogger("CountDL")
logger.setLevel(logging.DEBUG)

# create console handler with a higher log level
ch = logging.StreamHandler()
if args.debug:
    ch.setLevel(logging.DEBUG)
else:
    ch.setLevel(logging.INFO)
# create formatter and add it to the handlers
formatter_ch = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter_ch)
# add the handlers to the logger
logger.addHandler(ch)
logger.debug(sys.argv)




sptree_fn=args.sp_tree

rec_trees =glob.glob("%s/*" % args.rec_trees_dir)

Res_dict = {}
TotalD = 0

sp_t = Tree(sptree_fn,format=1)
i=0
for n in sp_t.traverse():
    Res_dict[str(i)] = {"Sp": n.name, "D":0, "S":0}
    i+=1

for tree_fn in rec_trees:
    print tree_fn
    t = Tree(tree_fn.strip())
    for n in t.traverse():
        Res_dict[str(n.S)][n.Ev] +=1

for k in Res_dict.keys():
    TotalD += Res_dict[k]["D"]

with open(args.o+".nbD.txt", "w") as f:
    f.write(str(TotalD) + "\n")

i=0
for n in sp_t.traverse():
    n.add_features(D=Res_dict[str(i)]["D"])
    n.add_features(S=Res_dict[str(i)]["S"])
    i+=1
    
sp_t.write(format=1, features=["D", "S"],outfile=args.o+".nhx", format_root_node=True)


from ete3 import TreeStyle, TextFace
# Basic tree style
tree_style = TreeStyle()
tree_style.show_leaf_name = False
tree_style.show_branch_length = False
tree_style.show_scale = False
tree_style.extra_branch_line_type=0 # 0=solid, 1=dashed, 2=dotted
tree_style.extra_branch_line_color="black"

i=0
for n in sp_t.traverse():
    n.dist=5
    Dstrinq = " D: " + str(Res_dict[str(i)]["D"])
    Sstrinq = " S: " + str(Res_dict[str(i)]["S"])
    n.add_face(TextFace(Dstrinq, fgcolor="#1C01AC"), column=0, position = 'branch-bottom')
    n.add_face(TextFace(Sstrinq, fgcolor="#800080"), column=0, position = 'branch-bottom')
    if n.is_leaf():
        n.add_face(TextFace(" "+n.name, fsize=12, fstyle="italic",), column=0,  position = "branch-right")
    elif n.name:
        n.add_face(TextFace(" "+n.name, fsize=12, fstyle="italic",), column=0,  position = "branch-top")
    i+=1

sp_t.render(args.o+".svg", tree_style=tree_style)


