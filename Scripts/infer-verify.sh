# get back to root
cd ..

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

cd ..

# compile and run hypertraps code
gcc -o3 hypertraps.c -lm -o hypertraps.ce

# two-pathway ("cross") case studies
./hypertraps.ce --obs VerifyData/synth-cross-samples-0.txt --starttimes VerifyData/synth-cross-times-0.txt --seed 1 --length 4 --kernel 5 --label VerifyData/cross-0 > VerifyData/cross03125.tmp &
./hypertraps.ce --obs VerifyData/synth-cross-samples-1.txt --starttimes VerifyData/synth-cross-times-1.txt --seed 1 --length 4 --kernel 5 --label VerifyData/cross-1 > VerifyData/cross13125.tmp &
./hypertraps.ce --obs VerifyData/synth-cross-samples-2.txt --starttimes VerifyData/synth-cross-times-2.txt --seed 1 --length 4 --kernel 5 --label VerifyData/cross-2 > VerifyData/cross23125.tmp &

# specific L=3 cubes with easy and hard parameterisations
./hypertraps.ce --obs VerifyData/synth-easycube-data.txt --starttimes VerifyData/synth-easycube-time.txt --seed 1 --length 4 --kernel 5 --label VerifyData/easycube > VerifyData/hardcube3125.tmp &
./hypertraps.ce --obs VerifyData/synth-hardcube-data.txt --starttimes VerifyData/synth-hardcube-time.txt --seed 1 --length 4 --kernel 5 --label VerifyData/hardcube > VerifyData/hardcube3125.tmp &

