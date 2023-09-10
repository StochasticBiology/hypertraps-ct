mkdir VerifyData/
cd Verify
gcc -o3 generate-cross.c -lm -o generate-cross.ce
./generate-cross.ce
cp synth* ../VerifyData

cd ../Inference

gcc -o3 hypertraps-all.c -lm -o hypertraps-all.ce

./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-0.txt --length 4 --label ../VerifyData/test-cross-1 > ../VerifyData/t-simple-1.tmp &
./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-0.txt --length 4 --label ../VerifyData/test-cross-2 --walkers 2000 > ../VerifyData/t-simple-2.tmp &
./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-0.txt --length 4 --label ../VerifyData/test-cross-3 --walkers 20 > ../VerifyData/t-simple-3.tmp &
./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-0.txt --length 4 --label ../VerifyData/test-cross-sa --sa > ../VerifyData/t-simple-sa.tmp &
./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-0.txt --length 4 --label ../VerifyData/test-cross-sgd --sgd > ../VerifyData/t-simple-sgd.tmp &
./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-0.txt --times ../VerifyData/synth-cross-times-0.txt --length 4 --label ../VerifyData/test-cross-ct-1 > ../VerifyData/t-simple-ct-1.tmp &
./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-0.txt --times ../VerifyData/synth-cross-times-0.txt --length 4 --walkers 2000 --label ../VerifyData/test-cross-ct-2 > ../VerifyData/t-simple-ct-2.tmp &
./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-0.txt --times ../VerifyData/synth-cross-times-time.txt --length 4 --walkers 20 --label ../VerifyData/test-cross-ct-3 > ../VerifyData/t-simple-ct-3.tmp &
