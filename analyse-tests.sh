cd Inference

gcc -o3 posteriors.c -lm -o posteriors.ce

./posteriors.ce --posterior ../VerifyData/test-cross-1-posterior.txt > ../VerifyData/t-simple-1-a.tmp &
./posteriors.ce --posterior  ../VerifyData/test-cross-2-posterior.txt  > ../VerifyData/t-simple-2-a.tmp &
./posteriors.ce --posterior  ../VerifyData/test-cross-3-posterior.txt  > ../VerifyData/t-simple-3-a.tmp &
./posteriors.ce --posterior  ../VerifyData/test-cross-ct-1-posterior.txt > ../VerifyData/t-simple-ct-1-a.tmp &
./posteriors.ce --posterior  ../VerifyData/test-cross-ct-2-posterior.txt > ../VerifyData/t-simple-ct-2-a.tmp &
./posteriors.ce --posterior  ../VerifyData/test-cross-ct-3-posterior.txt > ../VerifyData/t-simple-ct-3-a.tmp &
