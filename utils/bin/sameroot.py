#!/usr/bin/python
# coding: utf-8

#  sameroot.py
#
#  Copyright 2017 Carine Rey <carine.rey@ens-lyon.fr>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#

import sys
from ete3 import Tree

trees_fn = sys.argv[1]
trees_sameroot_fn = sys.argv[2]

F = open(trees_fn, "r")

tree_l = F.read().split("\n")
F.close()

root = None
trees_sameroot_s = []
for tree in tree_l:
    if tree:
        print tree
        t =  Tree(tree)
        if not root:
          root = t.get_leaf_names()[0]
          t.set_outgroup(root)
        else:
          t.set_outgroup(root)
        trees_sameroot_s.append(t.write(format=9))

F = open(trees_sameroot_fn, "w")
F.write("\n".join(trees_sameroot_s) + "\n")
F.close()
