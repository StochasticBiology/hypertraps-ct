# perform various more technical tests about optimisation and convergence for different sampling and optimisation schemes

mkdir VerifyData/
cd Verify
gcc -o3 generate-cross.c -lm -o generate-cross.ce
./generate-cross.ce
gcc -o3 generate-xor.c -lm -o generate-xor.ce
./generate-xor.ce

cp synth* ../VerifyData

cd ../Inference

gcc -o3 hypertraps-all.c -lm -o hypertraps-all.ce

## discrete time
# different sampling walkers
./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-0.txt --length 4 --label ../VerifyData/test-cross-1 > ../VerifyData/t-simple-1.tmp &
./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-0.txt --length 4 --label ../VerifyData/test-cross-2 --walkers 2000 > ../VerifyData/t-simple-2.tmp &
./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-0.txt --length 4 --label ../VerifyData/test-cross-3 --walkers 20 > ../VerifyData/t-simple-3.tmp &
# different optimisers
./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-0.txt --length 4 --label ../VerifyData/test-cross-sa --sa > ../VerifyData/t-simple-sa.tmp &
./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-0.txt --length 4 --label ../VerifyData/test-cross-sgd --sgd > ../VerifyData/t-simple-sgd.tmp &
# different parameter structures
./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-0.txt --length 4 --label ../VerifyData/test-cross-mod--1 --model -1 > ../VerifyData/t-simple-mod--1.tmp &
./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-0.txt --length 4 --label ../VerifyData/test-cross-mod-0 --model 0 > ../VerifyData/t-simple-mod-0.tmp &
./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-0.txt --length 4 --label ../VerifyData/test-cross-mod-1 --model 1 > ../VerifyData/t-simple-mod-1.tmp &
./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-0.txt --length 4 --label ../VerifyData/test-cross-mod-2 --model 2 > ../VerifyData/t-simple-mod-2.tmp &
./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-0.txt --length 4 --label ../VerifyData/test-cross-mod-3 --model 3 > ../VerifyData/t-simple-mod-3.tmp &
./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-0.txt --length 4 --label ../VerifyData/test-cross-mod-4 --model 4 > ../VerifyData/t-simple-mod-4.tmp &

### discrete time XOR problem
./hypertraps-all.ce --obs ../VerifyData/synth-xor.txt --length 4 --label ../VerifyData/test-xor-mod--1 --model -1 > ../VerifyData/t-xor-mod--1.tmp &
./hypertraps-all.ce --obs ../VerifyData/synth-xor.txt --length 4 --label ../VerifyData/test-xor-mod-0 --model 0 > ../VerifyData/t-xor-mod-0.tmp &
./hypertraps-all.ce --obs ../VerifyData/synth-xor.txt --length 4 --label ../VerifyData/test-xor-mod-1 --model 1 > ../VerifyData/t-xor-mod-1.tmp &
./hypertraps-all.ce --obs ../VerifyData/synth-xor.txt --length 4 --label ../VerifyData/test-xor-mod-2 --model 2 > ../VerifyData/t-xor-mod-2.tmp &
./hypertraps-all.ce --obs ../VerifyData/synth-xor.txt --length 4 --label ../VerifyData/test-xor-mod-3 --model 3 > ../VerifyData/t-xor-mod-3.tmp &
./hypertraps-all.ce --obs ../VerifyData/synth-xor.txt --length 4 --label ../VerifyData/test-xor-mod-4 --model 4 > ../VerifyData/t-xor-mod-4.tmp &


## continuous time
# different sampling walkers
./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-0.txt --times ../VerifyData/synth-cross-times-0.txt --length 4 --label ../VerifyData/test-cross-ct-1 > ../VerifyData/t-simple-ct-1.tmp &
./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-0.txt --times ../VerifyData/synth-cross-times-0.txt --length 4 --walkers 2000 --label ../VerifyData/test-cross-ct-2 > ../VerifyData/t-simple-ct-2.tmp &
./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-0.txt --times ../VerifyData/synth-cross-times-0.txt --length 4 --walkers 20 --label ../VerifyData/test-cross-ct-3 > ../VerifyData/t-simple-ct-3.tmp &
# different optimisers
./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-0.txt --times ../VerifyData/synth-cross-times-0.txt --length 4 --label ../VerifyData/test-cross-ct-sa --sa > ../VerifyData/t-simple-ct-sa.tmp &
./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-0.txt --times ../VerifyData/synth-cross-times-0.txt --length 4 --label ../VerifyData/test-cross-ct-sgd --sgd > ../VerifyData/t-simple-ct-sgd.tmp &
# different parameter structures
./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-0.txt --times ../VerifyData/synth-cross-times-0.txt --length 4 --label ../VerifyData/test-cross-ct-mod--1 --model -1 > ../VerifyData/t-simple-ct-mod--1.tmp &
./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-0.txt --times ../VerifyData/synth-cross-times-0.txt --length 4 --label ../VerifyData/test-cross-ct-mod-0 --model 0 > ../VerifyData/t-simple-ct-mod-0.tmp &
./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-0.txt --times ../VerifyData/synth-cross-times-0.txt --length 4 --label ../VerifyData/test-cross-ct-mod-1 --model 1 > ../VerifyData/t-simple-ct-mod-1.tmp &
./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-0.txt --times ../VerifyData/synth-cross-times-0.txt --length 4 --label ../VerifyData/test-cross-ct-mod-2 --model 2 > ../VerifyData/t-simple-ct-mod-2.tmp &
./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-0.txt --times ../VerifyData/synth-cross-times-0.txt --length 4 --label ../VerifyData/test-cross-ct-mod-3 --model 3 > ../VerifyData/t-simple-ct-mod-3.tmp &
./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-0.txt --times ../VerifyData/synth-cross-times-0.txt --length 4 --label ../VerifyData/test-cross-ct-mod-4 --model 4 > ../VerifyData/t-simple-ct-mod-4.tmp &

