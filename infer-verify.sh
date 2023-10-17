# generate synthetic datasets
mkdir VerifyData/
cd Verify
# these various small bits of code generate synthetic case study data for verification
gcc -o3 generate.c -lm -o generate.ce
gcc -o3 generate-easycube.c -lm -o generate-easycube.ce
gcc -o3 generate-hardcube.c -lm -o generate-hardcube.ce
gcc -o3 generate-cross.c -lm -o generate-cross.ce
gcc -o3 randomcube.c -lm -o randomcube.ce
./generate.ce 78 > generated.tmp
./generate-easycube.ce
./generate-hardcube.ce
./generate-cross.ce
./randomcube.ce
cp synth* ../VerifyData

cd ../Inference

# compile and run hypertraps code
gcc -o3 hypertraps-all.c -lm -o hypertraps-all.ce

# these old case studies don't feature in the current writeup
#./hypertraps-all.ce --obs ../VerifyData/synth-0-data.txt --times ../VerifyData/synth-0-times.txt --seed 1 --length 5 --kernel 3 --label ../VerifyData/synth-0-1 > ../VerifyData/tmp-synth-1 &
#./hypertraps-all.ce --obs ../VerifyData/synth-0-data.txt --times ../VerifyData/synth-0-times.txt --seed 2 --length 5 --kernel 3 --label ../VerifyData/synth-0-2 > ../VerifyData/tmp-synth-2 &

#./hypertraps-all.ce --obs ../VerifyData/synth-1-data.txt --times ../VerifyData/synth-1-times.txt --seed 1 --length 5 --kernel 3 --label ../VerifyData/synth-1-1 > ../VerifyData/tmp-synth-3 &
#./hypertraps-all.ce --obs ../VerifyData/synth-1-data.txt --times ../VerifyData/synth-1-times.txt --seed 2 --length 5 --kernel 3 --label ../VerifyData/synth-1-2 > ../VerifyData/tmp-synth-4 &

#./hypertraps-all.ce --obs ../VerifyData/synth-2-data.txt --times ../VerifyData/synth-2-times.txt --seed 1 --length 5 --kernel 3 --label ../VerifyData/synth-2-1 > ../VerifyData/tmp-synth-5 &
#./hypertraps-all.ce --obs ../VerifyData/synth-2-data.txt --times ../VerifyData/synth-2-times.txt --seed 2 --length 5 --kernel 3 --label ../VerifyData/synth-2-2 > ../VerifyData/tmp-synth-6 &

# two-pathway ("cross") case studies
./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-0.txt --times ../VerifyData/synth-cross-times-0.txt --seed 1 --length 4 --kernel 5 --label ../VerifyData/cross-0 > ../VerifyData/cross03125.tmp &
./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-1.txt --times ../VerifyData/synth-cross-times-1.txt --seed 1 --length 4 --kernel 5 --label ../VerifyData/cross-1 > ../VerifyData/cross13125.tmp &
./hypertraps-all.ce --obs ../VerifyData/synth-cross-samples-2.txt --times ../VerifyData/synth-cross-times-2.txt --seed 1 --length 4 --kernel 5 --label ../VerifyData/cross-2 > ../VerifyData/cross23125.tmp &

# specific L=3 cubes with easy and hard parameterisations
./hypertraps-all.ce --obs ../VerifyData/synth-easycube-data.txt --times ../VerifyData/synth-easycube-time.txt --seed 1 --length 4 --kernel 5 --label ../VerifyData/easycube > ../VerifyData/hardcube3125.tmp &
./hypertraps-all.ce --obs ../VerifyData/synth-hardcube-data.txt --times ../VerifyData/synth-hardcube-time.txt --seed 1 --length 4 --kernel 5 --label ../VerifyData/hardcube > ../VerifyData/hardcube3125.tmp &

