cd Inference

gcc -o3 posteriors.c -lm -o posteriors.ce

## discrete time
# different sampling walkers
./posteriors.ce --posterior ../VerifyData/test-cross-1-posterior.txt > ../VerifyData/t-simple-1-a.tmp &
./posteriors.ce --posterior  ../VerifyData/test-cross-2-posterior.txt  > ../VerifyData/t-simple-2-a.tmp &
./posteriors.ce --posterior  ../VerifyData/test-cross-3-posterior.txt  > ../VerifyData/t-simple-3-a.tmp &
# different parameter structures
./posteriors.ce --posterior ../VerifyData/test-cross-mod--1-posterior.txt --model -1 > ../VerifyData/t-sample-mod--1-a.tmp &
./posteriors.ce --posterior ../VerifyData/test-cross-mod-0-posterior.txt --model 0 > ../VerifyData/t-sample-mod-0-a.tmp &
./posteriors.ce --posterior ../VerifyData/test-cross-mod-1-posterior.txt --model 1 > ../VerifyData/t-sample-mod-1-a.tmp &
./posteriors.ce --posterior ../VerifyData/test-cross-mod-2-posterior.txt --model 2 > ../VerifyData/t-sample-mod-2-a.tmp &
./posteriors.ce --posterior ../VerifyData/test-cross-mod-3-posterior.txt --model 3 > ../VerifyData/t-sample-mod-3-a.tmp &
./posteriors.ce --posterior ../VerifyData/test-cross-mod-4-posterior.txt --model 4 > ../VerifyData/t-sample-mod-4-a.tmp &

### discrete time XOR problem
./posteriors.ce --posterior ../VerifyData/test-xor-mod--1-posterior.txt --model -1 > ../VerifyData/t-xor-mod--1-a.tmp &
./posteriors.ce --posterior ../VerifyData/test-xor-mod-0-posterior.txt --model 0 > ../VerifyData/t-xor-mod-0-a.tmp &
./posteriors.ce --posterior ../VerifyData/test-xor-mod-1-posterior.txt --model 1 > ../VerifyData/t-xor-mod-1-a.tmp &
./posteriors.ce --posterior ../VerifyData/test-xor-mod-2-posterior.txt --model 2 > ../VerifyData/t-xor-mod-2-a.tmp &
./posteriors.ce --posterior ../VerifyData/test-xor-mod-3-posterior.txt --model 3 > ../VerifyData/t-xor-mod-3-a.tmp &
./posteriors.ce --posterior ../VerifyData/test-xor-mod-4-posterior.txt --model 4 > ../VerifyData/t-xor-mod-4-a.tmp &

## continuous time
# different sampling walkers
./posteriors.ce --posterior  ../VerifyData/test-cross-ct-1-posterior.txt > ../VerifyData/t-simple-ct-1-a.tmp &
./posteriors.ce --posterior  ../VerifyData/test-cross-ct-2-posterior.txt > ../VerifyData/t-simple-ct-2-a.tmp &
./posteriors.ce --posterior  ../VerifyData/test-cross-ct-3-posterior.txt > ../VerifyData/t-simple-ct-3-a.tmp &
# different parameter structures
./posteriors.ce --posterior ../VerifyData/test-cross-ct-mod--1-posterior.txt --model -1 > ../VerifyData/t-sample-ct-mod--1-a.tmp &
./posteriors.ce --posterior ../VerifyData/test-cross-ct-mod-0-posterior.txt --model 0 > ../VerifyData/t-sample-ct-mod-0-a.tmp &
./posteriors.ce --posterior ../VerifyData/test-cross-ct-mod-1-posterior.txt --model 1 > ../VerifyData/t-sample-ct-mod-1-a.tmp &
./posteriors.ce --posterior ../VerifyData/test-cross-ct-mod-2-posterior.txt --model 2 > ../VerifyData/t-sample-ct-mod-2-a.tmp &
./posteriors.ce --posterior ../VerifyData/test-cross-ct-mod-3-posterior.txt --model 3 > ../VerifyData/t-sample-ct-mod-3-a.tmp &
./posteriors.ce --posterior ../VerifyData/test-cross-ct-mod-4-posterior.txt --model 4 > ../VerifyData/t-sample-ct-mod-4-a.tmp &
