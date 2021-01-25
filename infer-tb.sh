cd Inference
gcc -o3 hypertraps-all.c -lm -o hypertraps-all.ce

# discrete time approaches
./hypertraps-all.ce ../Data/ng.2878-S2.txt-cooked.txt-data.txt 0 1 3 4 0 0 > ../Data/tb-dt-1.tmp &
./hypertraps-all.ce ../Data/ng.2878-S2.txt-cooked.txt-data.txt 0 2 3 4 0 0 > ../Data/tb-dt-2.tmp &

# continuous time approaches
./hypertraps-all.ce ../Data/ng.2878-S2.txt-cooked.txt-data.txt ../Data/ng.2878-S2.txt-cooked.txt-datatime.txt ../Data/ng.2878-S2.txt-cooked.txt-datatime.txt 1 4 4 0 0 > ../Data/tb-ct-1.tmp &
./hypertraps-all.ce ../Data/ng.2878-S2.txt-cooked.txt-data.txt ../Data/ng.2878-S2.txt-cooked.txt-datatime.txt ../Data/ng.2878-S2.txt-cooked.txt-datatime.txt 2 4 4 0 0 > ../Data/tb-ct-2.tmp &
./hypertraps-all.ce ../Data/ng.2878-S2.txt-cooked.txt-data.txt ../Data/ng.2878-S2.txt-cooked.txt-datatime.txt ../Data/ng.2878-S2.txt-cooked.txt-datatime.txt 1 5 4 0 0 > ../Data/tb-ct-a-1.tmp &
./hypertraps-all.ce ../Data/ng.2878-S2.txt-cooked.txt-data.txt ../Data/ng.2878-S2.txt-cooked.txt-datatime.txt ../Data/ng.2878-S2.txt-cooked.txt-datatime.txt 2 5 4 0 0 > ../Data/tb-ct-a-2.tmp &
./hypertraps-all.ce ../Data/ng.2878-S2.txt-cooked.txt-data.txt ../Data/ng.2878-S2.txt-cooked.txt-datatime.txt ../Data/ng.2878-S2.txt-cooked.txt-datatime.txt 1 4 5 0 0 > ../Data/tb-ct-b-1.tmp &
./hypertraps-all.ce ../Data/ng.2878-S2.txt-cooked.txt-data.txt ../Data/ng.2878-S2.txt-cooked.txt-datatime.txt ../Data/ng.2878-S2.txt-cooked.txt-datatime.txt 2 4 5 0 0 > ../Data/tb-ct-b-2.tmp &
