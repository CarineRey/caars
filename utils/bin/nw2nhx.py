#!/usr/bin/python
# coding: utf-8

import sys
from ete2 import Tree

t = Tree(sys.argv[1])
print t.write(format=1, features=["support"])
