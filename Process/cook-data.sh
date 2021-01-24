# needs two arguments: a tree and a set of barcodes
# the tree should be in Newick format
# the barcodes should be in CSV format with the sample name as the first column, with one header row

# clean up the tree for use with software
tr -d '\n' < $1 | sed 's/ /_/g' | sed "s/'//g" > $1-cleaned.phy

# prune the tree to just include the nodes in the barcodes file
python prune-tree.py $1-cleaned.phy $2

# introduce dummy internal labels to the pruned tree
gcc internal-labels.c -o internal-labels.ce
./internal-labels.ce $1-clean.phy-pruned.txt > $1-cooked.txt

# construct ancestral states and hence transition datafiles
python parse-new.py $1-cooked.txt $2


