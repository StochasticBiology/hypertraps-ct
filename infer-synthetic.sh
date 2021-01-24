# generate synthetic datasets
cd Tests
gcc -o3 generate.c -lm -o generate.ce
./generate.ce 78 > generated.tmp

# copy these to inference folder, and follow them
cp synth* ../Inference
cd ../Inference

# compile and run hypertraps code
gcc -o3 hypertraps-all.c -lm -o hypertraps-all.ce
./hypertraps-all.ce synth-0-data.txt synth-0-times.txt synth-0-times.txt 1 3 3 0 0 > tmp-synth-1 &
./hypertraps-all.ce synth-0-data.txt synth-0-times.txt synth-0-times.txt 2 3 3 0 0 > tmp-synth-2 &

./hypertraps-all.ce synth-1-data.txt synth-1-times.txt synth-1-times.txt 1 3 3 0 0 > tmp-synth-3 &
./hypertraps-all.ce synth-1-data.txt synth-1-times.txt synth-1-times.txt 2 3 3 0 0 > tmp-synth-4 &

./hypertraps-all.ce synth-2-data.txt synth-2-times.txt synth-2-times.txt 1 3 3 0 0 > tmp-synth-5 &
./hypertraps-all.ce synth-2-data.txt synth-2-times.txt synth-2-times.txt 2 3 3 0 0 > tmp-synth-6 &
