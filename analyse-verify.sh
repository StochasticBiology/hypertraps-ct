# get to inference folder
cd Inference

# compile and run posterior analysis code
gcc -o3 posteriors.c -lm -o posteriors.ce
#./posteriors.ce --posterior synth-0-1-posterior.txt > ../VerifyData/posterior-synth-0-1.tmp &
#./posteriors.ce --posterior synth-0-2-posterior.txt > ../VerifyData/posterior-synth-0-2.tmp &
#./posteriors.ce --posterior synth-1-1-posterior.txt > ../VerifyData/posterior-synth-1-1.tmp &
#./posteriors.ce --posterior synth-1-2-posterior.txt > ../VerifyData/posterior-synth-1-2.tmp &
#./posteriors.ce --posterior synth-2-1-posterior.txt > ../VerifyData/posterior-synth-2-1.tmp &
#./posteriors.ce --posterior synth-2-2-posterior.txt > ../VerifyData/posterior-synth-2-2.tmp &

./posteriors.ce --posterior ../VerifyData/easycube-posterior.txt > ../VerifyData/posterior-easycube.tmp &
./posteriors.ce --posterior ../VerifyData/hardcube-posterior.txt > ../VerifyData/posterior-hardcube.tmp &
./posteriors.ce --posterior ../VerifyData/cross-0-posterior.txt > ../VerifyData/posterior-cross-0.tmp &
./posteriors.ce --posterior ../VerifyData/cross-1-posterior.txt > ../VerifyData/posterior-cross-1.tmp &
./posteriors.ce --posterior ../VerifyData/cross-2-posterior.txt > ../VerifyData/posterior-cross-2.tmp &
