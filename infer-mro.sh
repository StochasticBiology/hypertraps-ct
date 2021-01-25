cd Inference
gcc -o3 hypertraps-all.c -lm -o hypertraps-all.ce

### NCBI, discrete time
./hypertraps-all.ce ../Data/mro-ncbi-tree.phy-cooked.txt-data.txt 0 1 3 7 1 0 > ../Data/mro-expt-1.txt &
./hypertraps-all.ce ../Data/mro-ncbi-tree.phy-cooked.txt-data.txt 0 2 3 7 1 0 > ../Data/mro-expt-2.txt &

### TimeTree, discrete time
./hypertraps-all.ce ../Data/mro-tree-tt-format.phy-data.txt 0 1 3 7 1 0 > ../Data/mro-expt-3.txt &
./hypertraps-all.ce ../Data/mro-tree-tt-format.phy-data.txt 0 2 3 7 1 0 > ../Data/mro-expt-4.txt &

### TimeTree, continuous time
./hypertraps-all.ce ../Data/mro-tree-tt-format.phy-data.txt ../Data/mro-tree-tt-format.phy-datatime.txt ../Data/mro-tree-tt-format.phy-datatime.txt 1 4 7 1 0 > ../Data/mro-expt-5.txt &
./hypertraps-all.ce ../Data/mro-tree-tt-format.phy-data.txt ../Data/mro-tree-tt-format.phy-datatime.txt ../Data/mro-tree-tt-format.phy-datatime.txt 2 4 7 1 0 > ../Data/mro-expt-6.txt &

### TimeTree+, discrete time
./hypertraps-all.ce ../Data/mro-tree-ttplus-format.phy-data.txt 0 1 3 7 1 0 > ../Data/mro-expt-7.txt &
./hypertraps-all.ce ../Data/mro-tree-ttplus-format.phy-data.txt 0 2 3 7 1 0 > ../Data/mro-expt-8.txt &

### TimeTree+, continuous time
./hypertraps-all.ce ../Data/mro-tree-ttplus-format.phy-data.txt ../Data/mro-ttplus-1.txt ../Data/mro-ttplus-2.txt 1 4 7 1 0 > ../Data/mro-expt-9.txt &
./hypertraps-all.ce ../Data/mro-tree-ttplus-format.phy-data.txt ../Data/mro-ttplus-1.txt ../Data/mro-ttplus-2.txt 2 4 7 1 0 > ../Data/mro-expt-10.txt &
