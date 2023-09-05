# get to inference folder
cd Inference

# compile and run posterior analysis code
gcc -o3 posteriors.c -lm -o posteriors.ce
./posteriors.ce 0 10 ../VerifyData/synth-0-data.txt-posterior-1-1-3-3-0.txt > ../VerifyData/posterior-synth-0-1.tmp &
./posteriors.ce 0 10 ../VerifyData/synth-0-data.txt-posterior-1-2-3-3-0.txt > ../VerifyData/posterior-synth-0-2.tmp &
./posteriors.ce 0 10 ../VerifyData/synth-1-data.txt-posterior-1-1-3-3-0.txt > ../VerifyData/posterior-synth-1-1.tmp &
./posteriors.ce 0 10 ../VerifyData/synth-1-data.txt-posterior-1-2-3-3-0.txt > ../VerifyData/posterior-synth-1-2.tmp &
./posteriors.ce 0 10 ../VerifyData/synth-2-data.txt-posterior-1-1-3-3-0.txt > ../VerifyData/posterior-synth-2-1.tmp &
./posteriors.ce 0 10 ../VerifyData/synth-2-data.txt-posterior-1-2-3-3-0.txt > ../VerifyData/posterior-synth-2-2.tmp &

./posteriors.ce 0 10 ../VerifyData/synth-easycube-data.txt-posterior-1-1-2-5-0.txt > ../VerifyData/posterior-easycube.tmp &
./posteriors.ce 0 10 ../VerifyData/synth-hardcube-data.txt-posterior-1-1-2-5-0.txt > ../VerifyData/posterior-hardcube.tmp &
./posteriors.ce 0 10 ../VerifyData/synth-cross-samples-0.txt-posterior-1-1-2-5-0.txt > ../VerifyData/posterior-cross-0.tmp &
./posteriors.ce 0 10 ../VerifyData/synth-cross-samples-1.txt-posterior-1-1-2-5-0.txt > ../VerifyData/posterior-cross-1.tmp &
./posteriors.ce 0 10 ../VerifyData/synth-cross-samples-2.txt-posterior-1-1-2-5-0.txt > ../VerifyData/posterior-cross-2.tmp &
