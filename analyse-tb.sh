cd Inference
gcc -o3 posteriors.c -lm -o posteriors.ce

# direct time
./posteriors.ce --posterior ../Data/tb-dt-1-posterior.txt > ../Data/tb-dt-1-a.tmp &
./posteriors.ce --posterior ../Data/tb-dt-2-posterior.txt > ../Data/tb-dt-2-a.tmp &

# continuous time
./posteriors.ce --posterior ../Data/tb-ct-1-posterior.txt > ../Data/tb-ct-1-a.tmp &
./posteriors.ce --posterior ../Data/tb-ct-2-posterior.txt > ../Data/tb-ct-2-a.tmp &
./posteriors.ce --posterior ../Data/tb-ct-a-1-posterior.txt > ../Data/tb-ct-a-1-a.tmp &
./posteriors.ce --posterior ../Data/tb-ct-a-2-posterior.txt > ../Data/tb-ct-a-2-a.tmp &
./posteriors.ce --posterior ../Data/tb-ct-b-1-posterior.txt > ../Data/tb-ct-b-1-a.tmp &
./posteriors.ce --posterior ../Data/tb-ct-b-2-posterior.txt > ../Data/tb-ct-b-2-a.tmp &
