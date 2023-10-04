# perform various more technical tests about optimisation and convergence for different sampling and optimisation schemes

mkdir VerifyData/
cd Verify
gcc -o3 generate-cross.c -lm -o generate-cross.ce
./generate-cross.ce

cp synth* ../VerifyData
cp hi-order.txt ../VerifyData

cd ../Inference

gcc -o3 hypertraps-all.c -lm -o hypertraps-all.ce

## discrete time
# different sampling walkers
./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-1.txt --length 4 --label ../VerifyData/test-cross-1 > ../VerifyData/t-simple-1.tmp &
./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-1.txt --length 4 --label ../VerifyData/test-cross-2 --walkers 2000 > ../VerifyData/t-simple-2.tmp &
./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-1.txt --length 4 --label ../VerifyData/test-cross-3 --walkers 20 > ../VerifyData/t-simple-3.tmp &
# different optimisers
./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-1.txt --length 4 --label ../VerifyData/test-cross-sa --sa > ../VerifyData/t-simple-sa.tmp &
./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-1.txt --length 4 --label ../VerifyData/test-cross-sgd --sgd > ../VerifyData/t-simple-sgd.tmp &
# different parameter structures
./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-1.txt --length 4 --label ../VerifyData/test-cross-mod--1 --model -1 > ../VerifyData/t-simple-mod--1.tmp &
./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-1.txt --length 4 --label ../VerifyData/test-cross-mod-0 --model 0 > ../VerifyData/t-simple-mod-0.tmp &
./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-1.txt --length 4 --label ../VerifyData/test-cross-mod-1 --model 1 > ../VerifyData/t-simple-mod-1.tmp &
./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-1.txt --length 4 --label ../VerifyData/test-cross-mod-2 --model 2 > ../VerifyData/t-simple-mod-2.tmp &
./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-1.txt --length 4 --label ../VerifyData/test-cross-mod-3 --model 3 > ../VerifyData/t-simple-mod-3.tmp &
./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-1.txt --length 4 --label ../VerifyData/test-cross-mod-4 --model 4 > ../VerifyData/t-simple-mod-4.tmp &

### discrete time hi-order logic problem
./hypertraps-all.ce --obs ../VerifyData/hi-order.txt --length 4 --label ../VerifyData/test-ho-mod--1 --model -1 > ../VerifyData/t-ho-mod--1.tmp &
./hypertraps-all.ce --obs ../VerifyData/hi-order.txt --length 4 --label ../VerifyData/test-ho-mod-2 --model 2 > ../VerifyData/t-ho-mod-2.tmp &
./hypertraps-all.ce --obs ../VerifyData/hi-order.txt --length 4 --label ../VerifyData/test-ho-mod-3 --model 3 > ../VerifyData/t-ho-mod-3.tmp &
./hypertraps-all.ce --obs ../VerifyData/hi-order.txt --length 5 --label ../VerifyData/test-ho-mod--1-l --model -1 > ../VerifyData/t-ho-mod--1-l.tmp &
./hypertraps-all.ce --obs ../VerifyData/hi-order.txt --length 5 --label ../VerifyData/test-ho-mod-2-l --model 2 > ../VerifyData/t-ho-mod-2-l.tmp &
./hypertraps-all.ce --obs ../VerifyData/hi-order.txt --length 5 --label ../VerifyData/test-ho-mod-3-l --model 3 > ../VerifyData/t-ho-mod-3-l.tmp &
./hypertraps-all.ce --obs ../VerifyData/hi-order.txt --length 4 --sa --label ../VerifyData/test-ho-mod--1-sa --model -1 > ../VerifyData/t-ho-mod--1-sa.tmp &
./hypertraps-all.ce --obs ../VerifyData/hi-order.txt --length 4 --sa --label ../VerifyData/test-ho-mod-2-sa --model 2 > ../VerifyData/t-ho-mod-2-sa.tmp &
./hypertraps-all.ce --obs ../VerifyData/hi-order.txt --length 4 --sa --label ../VerifyData/test-ho-mod-3-sa --model 3 > ../VerifyData/t-ho-mod-3-sa.tmp &

## continuous time
# different sampling walkers
./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-1.txt --times ../VerifyData/synth-cross-times-1.txt --length 4 --label ../VerifyData/test-cross-ct-1 > ../VerifyData/t-simple-ct-1.tmp &
./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-1.txt --times ../VerifyData/synth-cross-times-1.txt --length 4 --walkers 2000 --label ../VerifyData/test-cross-ct-2 > ../VerifyData/t-simple-ct-2.tmp &
./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-1.txt --times ../VerifyData/synth-cross-times-1.txt --length 4 --walkers 20 --label ../VerifyData/test-cross-ct-3 > ../VerifyData/t-simple-ct-3.tmp &
# different optimisers
./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-1.txt --times ../VerifyData/synth-cross-times-1.txt --length 4 --label ../VerifyData/test-cross-ct-sa --sa > ../VerifyData/t-simple-ct-sa.tmp &
./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-1.txt --times ../VerifyData/synth-cross-times-1.txt --length 4 --label ../VerifyData/test-cross-ct-sgd --sgd > ../VerifyData/t-simple-ct-sgd.tmp &
# different parameter structures
./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-1.txt --times ../VerifyData/synth-cross-times-1.txt --length 4 --label ../VerifyData/test-cross-ct-mod--1 --model -1 > ../VerifyData/t-simple-ct-mod--1.tmp &
./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-1.txt --times ../VerifyData/synth-cross-times-1.txt --length 4 --label ../VerifyData/test-cross-ct-mod-0 --model 0 > ../VerifyData/t-simple-ct-mod-0.tmp &
./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-1.txt --times ../VerifyData/synth-cross-times-1.txt --length 4 --label ../VerifyData/test-cross-ct-mod-1 --model 1 > ../VerifyData/t-simple-ct-mod-1.tmp &
./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-1.txt --times ../VerifyData/synth-cross-times-1.txt --length 4 --label ../VerifyData/test-cross-ct-mod-2 --model 2 > ../VerifyData/t-simple-ct-mod-2.tmp &
./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-1.txt --times ../VerifyData/synth-cross-times-1.txt --length 4 --label ../VerifyData/test-cross-ct-mod-3 --model 3 > ../VerifyData/t-simple-ct-mod-3.tmp &
./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-1.txt --times ../VerifyData/synth-cross-times-1.txt --length 4 --label ../VerifyData/test-cross-ct-mod-4 --model 4 > ../VerifyData/t-simple-ct-mod-4.tmp &

