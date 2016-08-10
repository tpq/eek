#!/bin/bash

DB=qtdb
cd "physionet.org/physiobank/database/$DB"

# Select sinus rhythm IDs only
IDS=(sel16265 sel16272 sel16273 sel16420 sel16483 sel16539 sel16773 sel16786 sel16795 sel17453)
len=${#IDS[@]}

mkdir -p /home/thom/Dropbox/R/projects/_lib-eek/eek/data-raw/qtdb

# Use for loop read all IDs
for (( i=0; i<$len; i++ ));
do
  sample=${IDS[$i]}
  echo $sample

  # Note, 2-channel ECGs
  rdsamp -r $sample -p > "/home/thom/Dropbox/R/projects/_lib-eek/eek/data-raw/qtdb/$sample-phys.txt"

  # Save first pass PQRST identification
  rdann -r $sample -a qt1 > "/home/thom/Dropbox/R/projects/_lib-eek/eek/data-raw/qtdb/$sample-qt1.txt"

done
