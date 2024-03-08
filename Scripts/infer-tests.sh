#!/usr/bin/env bash
# perform various more technical tests about optimisation and convergence for different sampling and optimisation schemes

# get back to root
cd ..

mkdir VerifyData/
cd Verify
gcc -o3 generate-cross.c -lm -o generate-cross.ce
./generate-cross.ce
gcc -o3 generate-big-crosses.c -lm -o generate-big-crosses.ce
./generate-big-crosses.ce

cp synth* ../VerifyData
cp hi-order.txt ../VerifyData

cd ..

gcc -o3 hypertraps.c -lm -o hypertraps.ce

# Clears all command line arguments passed to the scripts, required to save pids
set --

## discrete time
# different sampling walkers
./hypertraps.ce --obs VerifyData/synth-cross-samples-1.txt --transitionformat --length 4 --label VerifyData/test-cross-1 > VerifyData/t-simple-1.tmp &
set -- "$@" $! # Append the previous process id to args
./hypertraps.ce --obs VerifyData/synth-cross-samples-1.txt --transitionformat --length 4 --label VerifyData/test-cross-2 --walkers 2000 > VerifyData/t-simple-2.tmp &
set -- "$@" $!
./hypertraps.ce --obs VerifyData/synth-cross-samples-1.txt --transitionformat --length 4 --label VerifyData/test-cross-3 --walkers 20 > VerifyData/t-simple-3.tmp &
set -- "$@" $!
# different optimisers
./hypertraps.ce --obs VerifyData/synth-cross-samples-1.txt --transitionformat --length 4 --label VerifyData/test-cross-sa --sa > VerifyData/t-simple-sa.tmp &
set -- "$@" $!
./hypertraps.ce --obs VerifyData/synth-cross-samples-1.txt --transitionformat --length 4 --label VerifyData/test-cross-sgd --sgd > VerifyData/t-simple-sgd.tmp &
set -- "$@" $!
# different parameter structures
./hypertraps.ce --obs VerifyData/synth-cross-samples-1.txt --transitionformat --length 4 --label VerifyData/test-cross-mod--1 --model -1 > VerifyData/t-simple-mod--1.tmp &
set -- "$@" $!
./hypertraps.ce --obs VerifyData/synth-cross-samples-1.txt --transitionformat --length 4 --label VerifyData/test-cross-mod-0 --model 0 > VerifyData/t-simple-mod-0.tmp &
set -- "$@" $!
./hypertraps.ce --obs VerifyData/synth-cross-samples-1.txt --transitionformat --length 4 --label VerifyData/test-cross-mod-1 --model 1 > VerifyData/t-simple-mod-1.tmp &
set -- "$@" $!
./hypertraps.ce --obs VerifyData/synth-cross-samples-1.txt --transitionformat --length 4 --label VerifyData/test-cross-mod-2 --model 2 > VerifyData/t-simple-mod-2.tmp &
set -- "$@" $!
./hypertraps.ce --obs VerifyData/synth-cross-samples-1.txt --transitionformat --length 4 --label VerifyData/test-cross-mod-3 --model 3 > VerifyData/t-simple-mod-3.tmp &
set -- "$@" $!
./hypertraps.ce --obs VerifyData/synth-cross-samples-1.txt --transitionformat --length 4 --label VerifyData/test-cross-mod-4 --model 4 > VerifyData/t-simple-mod-4.tmp &
set -- "$@" $!

## different parameter space sizes
# cross-cube observations
./hypertraps.ce --obs VerifyData/synth-bigcross-10-hard-samples.txt --transitionformat --length 4 --outputtransitions 0 --kernel 3 --label VerifyData/test-bigcross-hard-10 > VerifyData/t-bigcross-hard-10.tmp &
set -- "$@" $!
./hypertraps.ce --obs VerifyData/synth-bigcross-30-hard-samples.txt --transitionformat --length 4 --outputtransitions 0 --kernel 3 --label VerifyData/test-bigcross-hard-30 > VerifyData/t-bigcross-hard-30.tmp &
set -- "$@" $!
./hypertraps.ce --obs VerifyData/synth-bigcross-50-hard-samples.txt --transitionformat --length 4 --outputtransitions 0 --kernel 3 --label VerifyData/test-bigcross-hard-50 > VerifyData/t-bigcross-hard-50.tmp &
set -- "$@" $!
./hypertraps.ce --obs VerifyData/synth-bigcross-70-hard-samples.txt --transitionformat --length 4 --outputtransitions 0 --kernel 3 --label VerifyData/test-bigcross-hard-70 > VerifyData/t-bigcross-hard-70.tmp &
set -- "$@" $!
./hypertraps.ce --obs VerifyData/synth-bigcross-90-hard-samples.txt --transitionformat --length 4 --outputtransitions 0 --kernel 3 --label VerifyData/test-bigcross-hard-90 > VerifyData/t-bigcross-hard-90.tmp &
set -- "$@" $!
# cross-cube with different single-walker optimisers
./hypertraps.ce --obs VerifyData/synth-bigcross-90-hard-samples.txt --transitionformat --length 5 --outputtransitions 0 --kernel 3 --label VerifyData/test-bigcross-hard-90-1 --walkers 1 > VerifyData/t-bigcross-hard-90-1.tmp &
set -- "$@" $!
./hypertraps.ce --obs VerifyData/synth-bigcross-90-hard-samples.txt --transitionformat --length 4 --outputtransitions 0 --kernel 3 --label VerifyData/test-bigcross-hard-90-1-sa --sa --walkers 1 > VerifyData/t-bigcross-hard-90-1-sa.tmp &
set -- "$@" $!
./hypertraps.ce --obs VerifyData/synth-bigcross-90-hard-samples.txt --transitionformat --length 4 --outputtransitions 0 --kernel 3 --label VerifyData/test-bigcross-hard-90-1-sgd --sgd --walkers 1 > VerifyData/t-bigcross-hard-90-1-sgd.tmp &
set -- "$@" $!

### regularisation demo
./hypertraps.ce --obs VerifyData/synth-cross-samples-1.txt -lscale 100 --sa --transitionformat --length 4 --label VerifyData/test-cross-mod--1-sa --model -1 --regularise > VerifyData/t-simple-mod--1-sa.tmp &
set -- "$@" $!
./hypertraps.ce --obs VerifyData/synth-cross-samples-1.txt -lscale 100 --sa --transitionformat --length 4 --label VerifyData/test-cross-mod-2-sa --model 2 --regularise > VerifyData/t-simple-mod-2-sa.tmp &
set -- "$@" $!
./hypertraps.ce --obs VerifyData/synth-cross-samples-1.txt -lscale 100 --sa --transitionformat --length 4 --label VerifyData/test-cross-mod-3-sa --model 3 --regularise > VerifyData/t-simple-mod-3-sa.tmp &
set -- "$@" $!

# PLI
./hypertraps.ce --obs VerifyData/synth-cross-samples-1.txt --transitionformat --length 4 --label VerifyData/test-cross-mod-2-pli --pli --walkers 200 > t-pli.tmp &
set -- "$@" $!

### discrete time hi-order logic problem
./hypertraps.ce --obs VerifyData/hi-order.txt --lscale 10 --transitionformat --length 4 --label VerifyData/test-ho-mod--1 --model -1 > VerifyData/t-ho-mod--1.tmp &
set -- "$@" $!
./hypertraps.ce --obs VerifyData/hi-order.txt --lscale 10 --transitionformat --length 4 --label VerifyData/test-ho-mod-2 --model 2 > VerifyData/t-ho-mod-2.tmp &
set -- "$@" $!
./hypertraps.ce --obs VerifyData/hi-order.txt --lscale 10 --transitionformat --length 4 --label VerifyData/test-ho-mod-3 --model 3 > VerifyData/t-ho-mod-3.tmp &
set -- "$@" $!
./hypertraps.ce --obs VerifyData/hi-order.txt --lscale 10 --transitionformat --length 5 --label VerifyData/test-ho-mod--1-l --model -1 > VerifyData/t-ho-mod--1-l.tmp &
set -- "$@" $!
./hypertraps.ce --obs VerifyData/hi-order.txt --lscale 10 --transitionformat --length 5 --label VerifyData/test-ho-mod-2-l --model 2 > VerifyData/t-ho-mod-2-l.tmp &
set -- "$@" $!
./hypertraps.ce --obs VerifyData/hi-order.txt --lscale 10 --transitionformat --length 5 --label VerifyData/test-ho-mod-3-l --model 3 > VerifyData/t-ho-mod-3-l.tmp &
set -- "$@" $!
./hypertraps.ce --obs VerifyData/hi-order.txt --lscale 10 --transitionformat --length 5 --sa --label VerifyData/test-ho-mod--1-sa --model -1 > VerifyData/t-ho-mod--1-sa.tmp &
set -- "$@" $!
./hypertraps.ce --obs VerifyData/hi-order.txt --lscale 10 --transitionformat --length 5 --sa --label VerifyData/test-ho-mod-2-sa --model 2 > VerifyData/t-ho-mod-2-sa.tmp &
set -- "$@" $!
./hypertraps.ce --obs VerifyData/hi-order.txt --lscale 10 --transitionformat --length 5 --sa --label VerifyData/test-ho-mod-3-sa --model 3 > VerifyData/t-ho-mod-3-sa.tmp &
set -- "$@" $!

## continuous time
# different sampling walkers
./hypertraps.ce --obs VerifyData/synth-cross-samples-1.txt --starttimes VerifyData/synth-cross-times-1.txt --transitionformat --length 4 --label VerifyData/test-cross-ct-1 > VerifyData/t-simple-ct-1.tmp &
set -- "$@" $!
./hypertraps.ce --obs VerifyData/synth-cross-samples-1.txt --starttimes VerifyData/synth-cross-times-1.txt --transitionformat --length 4 --walkers 2000 --label VerifyData/test-cross-ct-2 > VerifyData/t-simple-ct-2.tmp &
set -- "$@" $!
./hypertraps.ce --obs VerifyData/synth-cross-samples-1.txt --starttimes VerifyData/synth-cross-times-1.txt --transitionformat --length 4 --walkers 20 --label VerifyData/test-cross-ct-3 > VerifyData/t-simple-ct-3.tmp &
set -- "$@" $!
# different optimisers
./hypertraps.ce --obs VerifyData/synth-cross-samples-1.txt --starttimes VerifyData/synth-cross-times-1.txt --transitionformat --length 4 --label VerifyData/test-cross-ct-sa --sa > VerifyData/t-simple-ct-sa.tmp &
set -- "$@" $!
./hypertraps.ce --obs VerifyData/synth-cross-samples-1.txt --starttimes VerifyData/synth-cross-times-1.txt --transitionformat --length 4 --label VerifyData/test-cross-ct-sgd --sgd > VerifyData/t-simple-ct-sgd.tmp &
set -- "$@" $!
# different parameter structures
./hypertraps.ce --obs VerifyData/synth-cross-samples-1.txt --starttimes VerifyData/synth-cross-times-1.txt --transitionformat --length 4 --label VerifyData/test-cross-ct-mod--1 --model -1 > VerifyData/t-simple-ct-mod--1.tmp &
set -- "$@" $!
./hypertraps.ce --obs VerifyData/synth-cross-samples-1.txt --starttimes VerifyData/synth-cross-times-1.txt --transitionformat --length 4 --label VerifyData/test-cross-ct-mod-0 --model 0 > VerifyData/t-simple-ct-mod-0.tmp &
set -- "$@" $!
./hypertraps.ce --obs VerifyData/synth-cross-samples-1.txt --starttimes VerifyData/synth-cross-times-1.txt --transitionformat --length 4 --label VerifyData/test-cross-ct-mod-1 --model 1 > VerifyData/t-simple-ct-mod-1.tmp &
set -- "$@" $!
./hypertraps.ce --obs VerifyData/synth-cross-samples-1.txt --starttimes VerifyData/synth-cross-times-1.txt --transitionformat --length 4 --label VerifyData/test-cross-ct-mod-2 --model 2 > VerifyData/t-simple-ct-mod-2.tmp &
set -- "$@" $!
./hypertraps.ce --obs VerifyData/synth-cross-samples-1.txt --starttimes VerifyData/synth-cross-times-1.txt --transitionformat --length 4 --label VerifyData/test-cross-ct-mod-3 --model 3 > VerifyData/t-simple-ct-mod-3.tmp &
set -- "$@" $!
./hypertraps.ce --obs VerifyData/synth-cross-samples-1.txt --starttimes VerifyData/synth-cross-times-1.txt --transitionformat --length 4 --label VerifyData/test-cross-ct-mod-4 --model 4 > VerifyData/t-simple-ct-mod-4.tmp &
set -- "$@" $!

## different parameter space sizes
# cross-cube observations
./hypertraps.ce --obs VerifyData/synth-bigcross-10-hard-samples.txt --starttimes VerifyData/synth-bigcross-10-hard-times.txt --transitionformat --length 4 --outputtransitions 0 --kernel 3 --label VerifyData/test-bigcross-hard-ct-10 > VerifyData/t-bigcross-hard-ct-10.tmp &
set -- "$@" $!
./hypertraps.ce --obs VerifyData/synth-bigcross-30-hard-samples.txt --starttimes VerifyData/synth-bigcross-30-hard-times.txt --transitionformat --length 4 --outputtransitions 0 --kernel 3 --label VerifyData/test-bigcross-hard-ct-30 > VerifyData/t-bigcross-hard-ct-30.tmp &
set -- "$@" $!
./hypertraps.ce --obs VerifyData/synth-bigcross-50-hard-samples.txt --starttimes VerifyData/synth-bigcross-50-hard-times.txt --transitionformat --length 4 --outputtransitions 0 --kernel 3 --label VerifyData/test-bigcross-hard-ct-50 > VerifyData/t-bigcross-hard-ct-50.tmp &
set -- "$@" $!
./hypertraps.ce --obs VerifyData/synth-bigcross-70-hard-samples.txt --starttimes VerifyData/synth-bigcross-70-hard-times.txt --transitionformat --length 4 --outputtransitions 0 --kernel 3 --label VerifyData/test-bigcross-hard-ct-70 > VerifyData/t-bigcross-hard-ct-70.tmp &
set -- "$@" $!
./hypertraps.ce --obs VerifyData/synth-bigcross-90-hard-samples.txt --starttimes VerifyData/synth-bigcross-90-hard-times.txt --transitionformat --length 4 --outputtransitions 0 --kernel 3 --label VerifyData/test-bigcross-hard-ct-90 > VerifyData/t-bigcross-hard-ct-90.tmp &
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
