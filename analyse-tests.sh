cd Inference

gcc -o3 posteriors.c -lm -o posteriors.ce

## discrete time
# different sampling walkers
./posteriors.ce --posterior ../VerifyData/test-cross-1-posterior.txt > ../VerifyData/t-simple-1-a.tmp &
./posteriors.ce --posterior  ../VerifyData/test-cross-2-posterior.txt  > ../VerifyData/t-simple-2-a.tmp &
./posteriors.ce --posterior  ../VerifyData/test-cross-3-posterior.txt  > ../VerifyData/t-simple-3-a.tmp &
# different parameter structures
./posteriors.ce --posterior ../VerifyData/test-cross-mod--1-posterior.txt --model -1 --verbose > ../VerifyData/t-sample-mod--1-a.tmp &
./posteriors.ce --posterior ../VerifyData/test-cross-mod-0-posterior.txt --model 0 --verbose > ../VerifyData/t-sample-mod-0-a.tmp &
./posteriors.ce --posterior ../VerifyData/test-cross-mod-1-posterior.txt --model 1 --verbose > ../VerifyData/t-sample-mod-1-a.tmp &
./posteriors.ce --posterior ../VerifyData/test-cross-mod-2-posterior.txt --model 2 --verbose > ../VerifyData/t-sample-mod-2-a.tmp &
./posteriors.ce --posterior ../VerifyData/test-cross-mod-3-posterior.txt --model 3 --verbose > ../VerifyData/t-sample-mod-3-a.tmp &
./posteriors.ce --posterior ../VerifyData/test-cross-mod-4-posterior.txt --model 4 --verbose > ../VerifyData/t-sample-mod-4-a.tmp &

### discrete time hi-order logic problem
./posteriors.ce --posterior ../VerifyData/test-ho-mod--1-posterior.txt --model -1 --verbose > ../VerifyData/t-ho-mod--1-a.tmp &
./posteriors.ce --posterior ../VerifyData/test-ho-mod-2-posterior.txt --model 2 --verbose > ../VerifyData/t-ho-mod-2-a.tmp &
./posteriors.ce --posterior ../VerifyData/test-ho-mod-3-posterior.txt --model 3 --verbose > ../VerifyData/t-ho-mod-3-a.tmp &
./posteriors.ce --posterior ../VerifyData/test-ho-mod--1-l-posterior.txt --model -1 --verbose > ../VerifyData/t-ho-mod--1-l-a.tmp &
./posteriors.ce --posterior ../VerifyData/test-ho-mod-2-l-posterior.txt --model 2 --verbose > ../VerifyData/t-ho-mod-2-l-a.tmp &
./posteriors.ce --posterior ../VerifyData/test-ho-mod-3-l-posterior.txt --model 3 --verbose > ../VerifyData/t-ho-mod-3-l-a.tmp &
./posteriors.ce --posterior ../VerifyData/test-ho-mod--1-sa-best.txt --model -1 --sims 1000 --verbose > ../VerifyData/t-ho-mod--1-sa-a.tmp &
./posteriors.ce --posterior ../VerifyData/test-ho-mod-2-sa-best.txt --model 2 --sims 1000 --verbose > ../VerifyData/t-ho-mod-2-sa-a.tmp &
./posteriors.ce --posterior ../VerifyData/test-ho-mod-3-sa-best.txt --model 3 --sims 1000 --verbose > ../VerifyData/t-ho-mod-3-sa-a.tmp &
./posteriors.ce --posterior ../VerifyData/test-ho-mod--1-sa-regularised.txt --model -1 --sims 1000 --verbose > ../VerifyData/t-ho-mod--1-sa-r-a.tmp &
./posteriors.ce --posterior ../VerifyData/test-ho-mod-2-sa-regularised.txt --model 2 --sims 1000 --verbose > ../VerifyData/t-ho-mod-2-sa-r-a.tmp &
./posteriors.ce --posterior ../VerifyData/test-ho-mod-3-sa-regularised.txt --model 3 --sims 1000 --verbose > ../VerifyData/t-ho-mod-3-sa-r-a.tmp &

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
