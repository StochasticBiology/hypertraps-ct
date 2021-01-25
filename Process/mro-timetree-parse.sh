# the TimeTree database doesn't hold all our species of interest
# this code addresses this with two approaches -- first, using close relatives; and second, manually adding species with uncertain divergence times

# first, prune the tree to a pre-compiled list of species it does contain
python3 prune-tree.py ../Data/Eukaryota_species.nwk ../Data/mro-species-for-timetree.txt

# in some cases we can approximate by relabelling close relatives with our species of interest
cp ../Data/Eukaryota_species.nwk-pruned.txt ../Data/Eukaryota_species_edited.nwk
sed -i 's/Reclinomonas_americana/Stygiella_incarcerata/g' ../Data/Eukaryota_species_edited.nwk
sed -i 's/Euglypha_rotunda/Brevimastigomonas_motovehiculus/g' ../Data/Eukaryota_species_edited.nwk
sed -i 's/Piromyces_communis/Piromyces_sp\./g' ../Data/Eukaryota_species_edited.nwk
sed -i 's/Neocallimastix_frontalis/Neocallimastix_sp\./g' ../Data/Eukaryota_species_edited.nwk
sed -i 's/Spironucleus_barkhanus/Spironucleus/g' ../Data/Eukaryota_species_edited.nwk

# now add internal labels
gcc internal-labels.c -o internal-labels.ce
./internal-labels.ce ../Data/Eukaryota_species_edited.nwk > ../Data/mro-tree-tt-format.phy

# parsing this tree gives us a reduced, but precisely specified, dataset
python3 parse-new.py ../Data/mro-tree-tt-format.phy ../Data/mro-barcodes.csv 0.001 > ../Data/mro-tree-tt-format.phy-output.txt 1


# in other cases we need to manually add species
# we do this by creating new branches from the internal nodes that we know (from NCBI) correspond to the common ancestors of the species to be added
# but we don't have timings for these, so we put them in with small fake "label" timings, which we then search-and-replace in the resulting timings file with the corresponding timings range, which is from 0 to the branch time from the common ancestor

cp ../Data/mro-tree-tt-format.phy ../Data/mro-tree-ttplus-format.phy
# for example, we know from NCBI that Orpinomyces branches from the internal node that will have been labelled AA, so we insert Orpinomyces with a dummy time label of "1" in this internal position
sed -i 's/)AA/,Orpinomyces_sp\.:1)AA/g' ../Data/mro-tree-ttplus-format.phy
sed -i 's/)AC/,Rozella:2)AC/g' ../Data/mro-tree-ttplus-format.phy
sed -i 's/)AC/,Mitosporidium:2)AC/g' ../Data/mro-tree-ttplus-format.phy
sed -i 's/)AC/,Encephalitozoon_cuniculi:2)AC/g' ../Data/mro-tree-ttplus-format.phy
sed -i 's/)AH/,Acanthamoeba_castellanii:3)AH/g' ../Data/mro-tree-ttplus-format.phy
sed -i 's/)AJ/,Psalteriomonas_lanterna:4)AJ/g' ../Data/mro-tree-ttplus-format.phy
sed -i 's/)AL/,Paratrimastix_pyriformis:5)AL/g' ../Data/mro-tree-ttplus-format.phy
sed -i 's/)AL/,Dysnectes_brevis:5)AL/g' ../Data/mro-tree-ttplus-format.phy
sed -i 's/)AQ/,Cryptosporidium_muris:6)AQ/g' ../Data/mro-tree-ttplus-format.phy
sed -i 's/)AR/,Vitrella_brassicaformis:7)AR/g' ../Data/mro-tree-ttplus-format.phy
sed -i 's/)AR/,Chromera_velia:7)AR/g' ../Data/mro-tree-ttplus-format.phy
sed -i 's/)AS/,Mikrocytos_mackini:8)AS/g' ../Data/mro-tree-ttplus-format.phy
sed -i 's/)AS/,Cantina_marsupialis:8)AS/g' ../Data/mro-tree-ttplus-format.phy
sed -i 's/)AS/,Blastocystis_sp.:8)AS/g' ../Data/mro-tree-ttplus-format.phy
sed -i 's/)AU/,Pygsuia_biforma:9)AU/g' ../Data/mro-tree-ttplus-format.phy

# parsing this tree gives us an expanded dataset for which some timings are uncertain
python3 parse-new.py ../Data/mro-tree-ttplus-format.phy ../Data/mro-barcodes.csv 0.001 > ../Data/mro-tree-ttplus-format.phy-output.txt 1

# now populate the time ranges for these added entries
cp ../Data/mro-tree-ttplus-format.phy-datatime.txt ../Data/mro-ttplus-1.txt
cp ../Data/mro-tree-ttplus-format.phy-datatime.txt ../Data/mro-ttplus-2.txt

# lower bound time is zero for all added species
sed -i 's/0\.00[0-9]/0/g' ../Data/mro-ttplus-1.txt

# upper bound time is time since branching for all added species
sed -i 's/0\.001/0.03375/g' ../Data/mro-ttplus-2.txt
sed -i 's/0\.002/0.987452/g' ../Data/mro-ttplus-2.txt
sed -i 's/0\.003/1.48/g' ../Data/mro-ttplus-2.txt
sed -i 's/0\.004/0.966737/g' ../Data/mro-ttplus-2.txt
sed -i 's/0\.005/1.4957639/g' ../Data/mro-ttplus-2.txt
sed -i 's/0\.006/0.817175/g' ../Data/mro-ttplus-2.txt
sed -i 's/0\.007/1.29/g' ../Data/mro-ttplus-2.txt
sed -i 's/0\.008/1.768418/g' ../Data/mro-ttplus-2.txt
sed -i 's/0\.009/2.100699/g' ../Data/mro-ttplus-2.txt
