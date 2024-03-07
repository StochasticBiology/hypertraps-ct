#!/usr/bin/env bash
# prepare and wrangle data for TB case study

# get back to root
cd ..

# shift raw data to folder for manipulation
mkdir Data/
cp RawData/* Data/

cd Process
chmod +x *.sh

# prepare TB data
./cook-data.sh ../Data/ng.2878-S2.txt ../Data/tuberculosis-v5-header-19-29.csv 1000 0
