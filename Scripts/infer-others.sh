#!/usr/bin/env bash
# run HyperTraPS for various existing scientific case studies

# get back to root
cd ..

gcc -o3 hypertraps.c -lm -o hypertraps.ce

# Clears all command line arguments passed to the scripts, required to save pids
set --

# C4 via MCMC, PLI, SA
./hypertraps.ce --obs Data/c4-curated.csv --length 4 --kernel 4 --label Data/c4-ht > Data/c4-ht.tmp &
set -- "$@" $! # Append the previous process id to args
./hypertraps.ce --obs Data/c4-curated.csv --length 4 --kernel 4 --pli --label Data/c4-pli --walkers 5000 > Data/c4-pli.tmp &
set -- "$@" $!
./hypertraps.ce --obs Data/c4-curated.csv --length 4 --kernel 4 --sa --label Data/c4-sa > Data/c4-sa.tmp &
set -- "$@" $!

# ovarian cancer via MCMC, PLI, SA
./hypertraps.ce --obs Data/ovarian.txt --length 4 --kernel 4 --label Data/ovarian-ht > Data/ovarian-ht.tmp &
set -- "$@" $!
./hypertraps.ce --obs Data/ovarian.txt --length 4 --kernel 4 --pli --label Data/ovarian-pli > Data/ovarian-pli.tmp &
set -- "$@" $!
./hypertraps.ce --obs Data/ovarian.txt --length 4 --kernel 4 --sa --label Data/ovarian-sa > Data/ovarian-sa.tmp &
set -- "$@" $!

# other others
./hypertraps.ce --obs Data/jallow_dataset_binary_with2s.csv --outputtransitions 0 --walkers 2 --kernel 3 --label Data/malaria > Data/malaria.tmp &
set -- "$@" $!
./hypertraps.ce --obs Data/jallow_dataset_binary_with2s.csv --outputtransitions 0 --walkers 2 --kernel 3 --label Data/malaria-sa --sa > Data/malaria-sa.tmp &
set -- "$@" $!
./hypertraps.ce --obs Data/total-observations.txt-trans.txt --label Data/tools --transitionformat--length 4 --kernel 4 > Data/tools.tmp &
set -- "$@" $!

echo Running $# instances of HyperTraPS as subprocesses
echo Waiting for all subprocesses to finish
echo PIDS: $@

# Wait for all subprocess created by this script to return
for job in $@
do
	wait $job
	echo $job finished
done
