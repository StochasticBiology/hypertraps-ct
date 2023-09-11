cd Inference
gcc -o3 posteriors.c -lm -o posteriors.ce

# direct time
./posteriors.ce --posterior ../Data/mro-1-posterior.txt > ../Data/mro-1-a.tmp &
./posteriors.ce --posterior ../Data/mro-2-posterior.txt > ../Data/mro-2-a.tmp &

# continuous time
./posteriors.ce --posterior ../Data/mro-3-posterior.txt > ../Data/mro-3-a.tmp &
./posteriors.ce --posterior ../Data/mro-4-posterior.txt > ../Data/mro-4-a.tmp &
./posteriors.ce --posterior ../Data/mro-5-posterior.txt > ../Data/mro-5-a.tmp &
./posteriors.ce --posterior ../Data/mro-6-posterior.txt > ../Data/mro-6-a.tmp &
./posteriors.ce --posterior ../Data/mro-7-posterior.txt > ../Data/mro-7-a.tmp &
./posteriors.ce --posterior ../Data/mro-8-posterior.txt > ../Data/mro-8-a.tmp &
./posteriors.ce --posterior ../Data/mro-9-posterior.txt > ../Data/mro-9-a.tmp &
./posteriors.ce --posterior ../Data/mro-10-posterior.txt > ../Data/mro-10-a.tmp &
