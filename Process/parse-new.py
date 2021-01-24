# takes two arguments -- a (pruned) Newick tree (argument 1) and a list of "barcode" observations
# uses model assumptions to label internal nodes
# outputs transitions and times to new files

from ete3 import Tree
import sys

# open file 1
# read in Newick tree -- assumes that it's all on one line
fp = open(str(sys.argv[1]), "r")
nw = fp.readline()
fp.close()

# interprets this content as a phylogenetic tree
tree = Tree(nw, format=1) # newick subformat 1 to read internal node names

# open file 2
# populate a dictionary, labelled by species names, with barcodes of traits
mydict = {}
with open(str(sys.argv[2])) as f:
    for line in f:
       elements = line.rstrip("\n").split(",")
       key = elements[0]
       val = "".join(elements[1:])
       mydict[key] = val
       ntraits = len(val)

# go through our tree from leaves up
for node in tree.traverse("postorder"):
    # A verbose output
    if not node.is_leaf() and node.name not in mydict: # that is, for all internal (ancestral) nodes
        # ones vector
        ref = ['1']*ntraits
        # go through children and assign 1s to ancestor if a descendant has them
        for childnode in node.children:
            for i, c in enumerate(mydict[childnode.name]):
                if c == '0':
                    ref[i] = '0'
        mydict[node.name] = "".join(ref)

# output reconstructed nodes and timings
fp = open(str(sys.argv[1])+"-data.txt", "w")
fptime = open(str(sys.argv[1])+"-datatime.txt", "w")
for node in tree.traverse("postorder"):
    # A verbose output
    if not node.is_leaf(): # For all internal nodes
        for childnode in node.children:
            print mydict[node.name], "(", node.name, ") -> ", mydict[childnode.name], "(", childnode.name, ")", "=", childnode.dist*1000
            print >>fp, " ".join(str(x) for x in list(mydict[childnode.name]))
            print >>fp, " ".join(str(x) for x in list(mydict[node.name]))
            print >>fptime, childnode.dist*1000
fp.close()
fptime.close()
