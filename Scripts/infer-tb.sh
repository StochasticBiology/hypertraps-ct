#!/usr/bin/env bash
# infer TB evolutionary dynamics
# ensure prepare-all.sh has been run to wrangle data

# get back to root
cd ..

gcc -o3 hypertraps.c -lm -o hypertraps.ce

# Clears all command line arguments passed to the scripts, required to save pids
set --

# discrete time approaches
./hypertraps.ce --obs Data/ng.2878-S2.txt-cooked.txt-data.txt --seed 1 --transitionformat --length 5 --kernel 4 --label Data/tb-dt-1 > Data/tb-dt-1.tmp &
set -- "$@" $! # Append the previous process id to args
./hypertraps.ce --obs Data/ng.2878-S2.txt-cooked.txt-data.txt --seed 2 --transitionformat --length 5 --kernel 4 --label Data/tb-dt-2 > Data/tb-dt-2.tmp &
set -- "$@" $!
# non-MCMC
./hypertraps.ce --obs Data/ng.2878-S2.txt-cooked.txt-data.txt --seed 1 --transitionformat --length 5 --kernel 4 --label Data/tb-dt-sa --sa > Data/tb-dt-sa.tmp &
set -- "$@" $!
./hypertraps.ce --obs Data/ng.2878-S2.txt-cooked.txt-data.txt --seed 1 --transitionformat --length 5 --kernel 4 --label Data/tb-dt-sgd --sgd > Data/tb-dt-sgd.tmp &
set -- "$@" $!

# continuous time approaches
./hypertraps.ce --obs Data/ng.2878-S2.txt-cooked.txt-data.txt --starttimes Data/ng.2878-S2.txt-cooked.txt-datatime.txt --seed 1 --transitionformat --length 6 --kernel 4 --label Data/tb-ct-1 > Data/tb-ct-1.tmp &
set -- "$@" $!
./hypertraps.ce --obs Data/ng.2878-S2.txt-cooked.txt-data.txt --starttimes Data/ng.2878-S2.txt-cooked.txt-datatime.txt --seed 2 --transitionformat --length 6 --kernel 4 --label Data/tb-ct-2 > Data/tb-ct-2.tmp &
set -- "$@" $!
./hypertraps.ce --obs Data/ng.2878-S2.txt-cooked.txt-data.txt --starttimes Data/ng.2878-S2.txt-cooked.txt-datatime.txt --seed 1 --transitionformat --length 7 --kernel 4 --label Data/tb-ct-a-1 > Data/tb-ct-a-1.tmp &
set -- "$@" $!
./hypertraps.ce --obs Data/ng.2878-S2.txt-cooked.txt-data.txt --starttimes Data/ng.2878-S2.txt-cooked.txt-datatime.txt --seed 2 --transitionformat --length 7 --kernel 4 --label Data/tb-ct-a-2 > Data/tb-ct-a-2.tmp &
set -- "$@" $!
./hypertraps.ce --obs Data/ng.2878-S2.txt-cooked.txt-data.txt --starttimes Data/ng.2878-S2.txt-cooked.txt-datatime.txt --seed 1 --transitionformat --length 6 --kernel 5 --label Data/tb-ct-b-1 > Data/tb-ct-b-1.tmp &
set -- "$@" $!
./hypertraps.ce --obs Data/ng.2878-S2.txt-cooked.txt-data.txt --starttimes Data/ng.2878-S2.txt-cooked.txt-datatime.txt --seed 2 --transitionformat --length 6 --kernel 5 --label Data/tb-ct-b-2 > Data/tb-ct-b-2.tmp &
set -- "$@" $!
# non-MCMC
./hypertraps.ce --obs Data/ng.2878-S2.txt-cooked.txt-data.txt --starttimes Data/ng.2878-S2.txt-cooked.txt-datatime.txt --seed 1 --transitionformat --length 5 --kernel 5 --label Data/tb-ct-sa --sa > Data/tb-ct-sa.tmp &
set -- "$@" $!
./hypertraps.ce --obs Data/ng.2878-S2.txt-cooked.txt-data.txt --starttimes Data/ng.2878-S2.txt-cooked.txt-datatime.txt --seed 1 --transitionformat --length 5 --kernel 5 --label Data/tb-ct-sgd --sgd > Data/tb-ct-sgd.tmp &
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
