cd Inference
gcc -o3 hypertraps-all.c -lm -o hypertraps-all.ce

./hypertraps-all.ce --obs ../Data/c4-curated.csv --crosssectional --length 4 --kernel 4 --label ../Data/c4-ht > ../Data/c4-ht.tmp &
./hypertraps-all.ce --obs ../Data/c4-curated.csv --crosssectional --length 4 --kernel 4 --pli --label ../Data/c4-pli > ../Data/c4-pli.tmp &
./hypertraps-all.ce --obs ../Data/c4-curated.csv --crosssectional --length 4 --kernel 4 --sa --label ../Data/c4-sa > ../Data/c4-sa.tmp &

./hypertraps-all.ce --obs ../Data/ovarian.txt --crosssectional --length 4 --kernel 4 --label ../Data/ovarian-ht > ../Data/ovarian-ht.tmp &
./hypertraps-all.ce --obs ../Data/ovarian.txt --crosssectional --length 4 --kernel 4 --pli --label ../Data/ovarian-pli > ../Data/ovarian-pli.tmp &
./hypertraps-all.ce --obs ../Data/ovarian.txt --crosssectional --length 4 --kernel 4 --sa --label ../Data/ovarian-sa > ../Data/ovarian-sa.tmp &
