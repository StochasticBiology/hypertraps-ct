# takes four arguments -- a (pruned) Newick tree (argument 1) and a list of "barcode" observations, and a time scaling factor (argument 3), and whether to consider loss/gain of traits (1/0)
# uses model assumptions to label internal nodes
# outputs transitions and times to new files

from ete3 import Tree, TreeStyle, NodeStyle, TextFace
import sys

arg1 = sys.argv[1]
arg2 = sys.argv[2]
arg3 = sys.argv[3]
arg4 = sys.argv[4]

# open file 1
# read in Newick tree -- assumes that it's all on one line
fp = open(str(arg1), "r")
nw = fp.readline()
fp.close()

# interprets this content as a phylogenetic tree
tree = Tree(nw, format=1) # newick subformat 1 to read internal node names

# open file 2
# populate a dictionary, labelled by species names, with barcodes of traits
mydict = {}
with open(str(arg2)) as f:
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
        # reference (end state) vector
        if int(arg4) == 0:
            ref = ['1']*ntraits
        else:
            ref = ['0']*ntraits
        # go through children and assign 1s to ancestor if a descendant has them
        for childnode in node.children:
            for i, c in enumerate(mydict[childnode.name]):
                if int(arg4) == 0:
                    if c == '0':
                        ref[i] = '0'
                else:
                    if c == '1':
                        ref[i] = '1'
        mydict[node.name] = "".join(ref)

# output reconstructed nodes and timings
fp = open(str(arg1)+"-data.txt", "w")
fptime = open(str(arg1)+"-datatime.txt", "w")
for node in tree.traverse("postorder"):
    # A verbose output
    if not node.is_leaf(): # For all internal nodes
        for childnode in node.children:
            print(mydict[node.name], "(", node.name, ") -> ", mydict[childnode.name], "(", childnode.name, ")", "=", childnode.dist*float(arg3))
            print(" ".join(str(x) for x in list(mydict[node.name])), file=fp)
            print(" ".join(str(x) for x in list(mydict[childnode.name])), file=fp)
            print(childnode.dist*float(arg3), file=fptime)

fp.close()
fptime.close()

# label nodes with barcodes for checking output
ts = TreeStyle()
ts.show_leaf_name = True
for leaf in tree.iter_leaves():
  thisleafcontent = TextFace(" ".join(str(x) for x in list(mydict[leaf.name])))
  leaf.add_face(thisleafcontent, 0, "aligned")

# output check tree to file
fname = str(arg1)+"-check.png"
tree.render(str(fname), w=800, tree_style=ts)

