# needs four arguments: a tree and a set of barcodes, a time scaling factor, and loss/gain (1/0)
# the tree should be in Newick format
# the barcodes should be in CSV format with the sample name as the first column, with one header row

# clean up the tree for use with software
tr -d '\n' < $1 | sed 's/ /_/g' | sed "s/'//g" > $1-cleaned.phy

# prune the tree to just include the nodes in the barcodes file
python3 prune-tree.py $1-cleaned.phy $2

# introduce dummy internal labels to the pruned tree
gcc internal-labels.c -o internal-labels.ce
./internal-labels.ce $1-cleaned.phy-pruned.txt > $1-cooked.txt

# construct ancestral states and hence transition datafiles
python3 parse-new.py $1-cooked.txt $2 $3 $4 > $1-output.txt



