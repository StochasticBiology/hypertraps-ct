# takes two arguments: a tree in Newick format (argument 1) and a file containing a set of sample names in a rather specific format (CSV, first column is sample names, first row is header (ignored))
# prunes the tree so that only those sample names are represented, and outputs to a new file

from ete3 import Tree
import sys

# get the list of sample names from argument 2
nodes = []
with open(str(sys.argv[2])) as f:
   for line in f:
     elements = line.split(",")
     nodes.append(elements[0])

# pop the first (header) row
nodes.pop(0)

# read the Newick tree
fname = str(sys.argv[1])
fp = open(fname, "r")
nw = fp.readline()
fp.close()

tree = Tree(nw, quoted_node_names=True, format=1) # newick subformat 1 to read internal node names

# prune the tree
tree.prune(nodes, preserve_branch_length=True)

# output to new file
fp = open(fname+"-pruned.txt", "w")
fp.write(tree.write())
fp.close()


