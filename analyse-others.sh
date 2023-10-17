cd Inference

gcc -o3 posteriors.c -lm -o posteriors.ce

# C4 via MCMC, PLI, SA
./posteriors.ce --posterior ../Data/c4-ht-posterior.txt --featurenames ../Data/c4-trait-names.txt > ../Data/c4-ht-a.tmp &
./posteriors.ce --posterior ../Data/c4-pli-posterior.txt --featurenames ../Data/c4-trait-names.txt  > ../Data/c4-pli-a.tmp &
./posteriors.ce --posterior ../Data/c4-sa-best.txt --featurenames ../Data/c4-trait-names.txt  > ../Data/c4-sa-a.tmp &

# ovarian cancer via MCMC, PLI, SA
./posteriors.ce --posterior ../Data/ovarian-ht-posterior.txt --featurenames ../Data/ovarian-names.txt > ../Data/ovarian-ht-a.tmp &
./posteriors.ce --posterior ../Data/ovarian-pli-posterior.txt --featurenames ../Data/ovarian-names.txt > ../Data/ovarian-pli-a.tmp &
./posteriors.ce --posterior ../Data/ovarian-sa-best.txt --featurenames ../Data/ovarian-names.txt > ../Data/ovarian-sa-a.tmp &

# C4 via MCMC, PLI, SA -- max lik
./posteriors.ce --posterior ../Data/c4-ht-best.txt --featurenames ../Data/c4-trait-names.txt > ../Data/c4-ht-b-a.tmp &
./posteriors.ce --posterior ../Data/c4-pli-best.txt --featurenames ../Data/c4-trait-names.txt  > ../Data/c4-pli-b-a.tmp &
./posteriors.ce --posterior ../Data/c4-sa-best.txt --featurenames ../Data/c4-trait-names.txt  > ../Data/c4-sa-b-a.tmp &
