# infer mitochondrion-related organelle evolutionary dynamics
# ensure prepare-all.sh has been run to wrangle data

# get back to root
cd ..

gcc -o3 hypertraps.c -lm -o hypertraps.ce

### NCBI, discrete time
./hypertraps.ce --obs Data/mro-ncbi-tree.phy-cooked.txt-data.txt --seed 1 --transitionformat --length 5 --kernel 7 --losses --label Data/mro-1 > Data/mro-expt-1.txt &
./hypertraps.ce --obs Data/mro-ncbi-tree.phy-cooked.txt-data.txt --seed 2 --transitionformat --length 5 --kernel 7 --losses --label Data/mro-2 > Data/mro-expt-2.txt &

### TimeTree, discrete time
./hypertraps.ce --obs Data/mro-tree-tt-format.phy-data.txt --seed 1 --transitionformat --length 5 --kernel 7 --losses --label Data/mro-3 > Data/mro-expt-3.txt &
./hypertraps.ce --obs Data/mro-tree-tt-format.phy-data.txt --seed 2 --transitionformat --length 5 --kernel 7 --losses --label Data/mro-4 > Data/mro-expt-4.txt &

### TimeTree, continuous time
./hypertraps.ce --obs Data/mro-tree-tt-format.phy-data.txt --starttimes Data/mro-tree-tt-format.phy-datatime.txt --seed 1 --transitionformat --length 6 --kernel 7 --losses --label Data/mro-5 > Data/mro-expt-5.txt &
./hypertraps.ce --obs Data/mro-tree-tt-format.phy-data.txt --starttimes Data/mro-tree-tt-format.phy-datatime.txt --seed 2 --transitionformat --length 6 --kernel 7 --losses --label Data/mro-6 > Data/mro-expt-6.txt &

### TimeTree+, discrete time
./hypertraps.ce --obs Data/mro-tree-ttplus-format.phy-data.txt --seed 1 --transitionformat --length 5 --kernel 7 --losses --label Data/mro-7 > Data/mro-expt-7.txt &
./hypertraps.ce --obs Data/mro-tree-ttplus-format.phy-data.txt --seed 2 --transitionformat --length 5 --kernel 7 --losses --label Data/mro-8 > Data/mro-expt-8.txt &

### TimeTree+, continuous time
./hypertraps.ce --obs Data/mro-tree-ttplus-format.phy-data.txt --starttimes Data/mro-ttplus-1.txt --endtimes Data/mro-ttplus-2.txt --seed 1 --transitionformat --length 6 --kernel 7 --losses --label Data/mro-9 > Data/mro-expt-9.txt &
./hypertraps.ce --obs Data/mro-tree-ttplus-format.phy-data.txt --starttimes Data/mro-ttplus-1.txt --endtimes Data/mro-ttplus-2.txt --seed 2 --transitionformat --length 6 --kernel 7 --losses --label Data/mro-10 > Data/mro-expt-10.txt &

# non-MCMC -- just SGD here
### NCBI, discrete time
./hypertraps.ce --obs Data/mro-ncbi-tree.phy-cooked.txt-data.txt --seed 1 --transitionformat --length 4 --kernel 7 --losses --label Data/mro-1-sgd --sgd > Data/mro-expt-1-sgd.txt &
### TimeTree, discrete time
./hypertraps.ce --obs Data/mro-tree-tt-format.phy-data.txt --seed 1 --transitionformat --length 4 --kernel 7 --losses --label Data/mro-3-sgd --sgd > Data/mro-expt-3-sgd.txt &
### TimeTree, continuous time
./hypertraps.ce --obs Data/mro-tree-tt-format.phy-data.txt --starttimes Data/mro-tree-tt-format.phy-datatime.txt --seed 1 --transitionformat --length 4 --kernel 7 --losses --label Data/mro-5-sgd --sgd > Data/mro-expt-5-sgd.txt &
### TimeTree+, discrete time
./hypertraps.ce --obs Data/mro-tree-ttplus-format.phy-data.txt --seed 1 --transitionformat --length 4 --kernel 7 --losses --label Data/mro-7-sgd --sgd > Data/mro-expt-7-sgd.txt &
### TimeTree+, continuous time
./hypertraps.ce --obs Data/mro-tree-ttplus-format.phy-data.txt --starttimes Data/mro-ttplus-1.txt --endtimes Data/mro-ttplus-2.txt --seed 1 --transitionformat --length 4 --kernel 7 --losses --label Data/mro-9-sgd --sgd > Data/mro-expt-9-sgd.txt &

