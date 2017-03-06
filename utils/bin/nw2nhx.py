#!/usr/bin/python
# coding: utf-8

import sys
from ete2 import Tree

t = Tree(sys.argv[1])
t.write(outfile=sys.argv[2],format=1, features=[])

