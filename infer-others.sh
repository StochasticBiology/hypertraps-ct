cd Inference
gcc -o3 hypertraps-all.c -lm -o hypertraps-all.ce

# C4 via MCMC, PLI, SA
./hypertraps-all.ce --obs ../Data/c4-curated.csv --crosssectional --length 4 --kernel 4 --label ../Data/c4-ht > ../Data/c4-ht.tmp &
./hypertraps-all.ce --obs ../Data/c4-curated.csv --crosssectional --length 4 --kernel 4 --pli --label ../Data/c4-pli --walkers 5000 > ../Data/c4-pli.tmp &
./hypertraps-all.ce --obs ../Data/c4-curated.csv --crosssectional --length 4 --kernel 4 --sa --label ../Data/c4-sa > ../Data/c4-sa.tmp &

# ovarian cancer via MCMC, PLI, SA
./hypertraps-all.ce --obs ../Data/ovarian.txt --crosssectional --length 4 --kernel 4 --label ../Data/ovarian-ht > ../Data/ovarian-ht.tmp &
./hypertraps-all.ce --obs ../Data/ovarian.txt --crosssectional --length 4 --kernel 4 --pli --label ../Data/ovarian-pli > ../Data/ovarian-pli.tmp &
./hypertraps-all.ce --obs ../Data/ovarian.txt --crosssectional --length 4 --kernel 4 --sa --label ../Data/ovarian-sa > ../Data/ovarian-sa.tmp &

# other others
./hypertraps-all.ce --obs ../Data/jallow_dataset_binary_with2s.csv --crosssectional --outputtransitions 0 --walkers 2 --label malaria
./hypertraps-all.ce --obs ../Data/total-observations.txt-trans.txt --label tools
