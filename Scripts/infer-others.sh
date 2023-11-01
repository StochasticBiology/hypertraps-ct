# run HyperTraPS for various existing scientific case studies

# get back to root
cd ..

gcc -o3 hypertraps.c -lm -o hypertraps.ce

# C4 via MCMC, PLI, SA
./hypertraps.ce --obs Data/c4-curated.csv --length 4 --kernel 4 --label Data/c4-ht > Data/c4-ht.tmp &
./hypertraps.ce --obs Data/c4-curated.csv --length 4 --kernel 4 --pli --label Data/c4-pli --walkers 5000 > Data/c4-pli.tmp &
./hypertraps.ce --obs Data/c4-curated.csv --length 4 --kernel 4 --sa --label Data/c4-sa > Data/c4-sa.tmp &

# ovarian cancer via MCMC, PLI, SA
./hypertraps.ce --obs Data/ovarian.txt --length 4 --kernel 4 --label Data/ovarian-ht > Data/ovarian-ht.tmp &
./hypertraps.ce --obs Data/ovarian.txt --length 4 --kernel 4 --pli --label Data/ovarian-pli > Data/ovarian-pli.tmp &
./hypertraps.ce --obs Data/ovarian.txt --length 4 --kernel 4 --sa --label Data/ovarian-sa > Data/ovarian-sa.tmp &

# other others
./hypertraps.ce --obs Data/jallow_dataset_binary_with2s.csv --outputtransitions 0 --walkers 2 --kernel 3 --label Data/malaria > Data/malaria.tmp &
./hypertraps.ce --obs Data/jallow_dataset_binary_with2s.csv --outputtransitions 0 --walkers 2 --kernel 3 --label Data/malaria-sa --sa > Data/malaria-sa.tmp &
./hypertraps.ce --obs Data/total-observations.txt-trans.txt --label Data/tools --transitionformat--length 4 --kernel 4 > Data/tools.tmp &
