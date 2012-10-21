#!/bin/bash
#
#
### SGE OPTIONS ###
#
#$ -cwd
#$ -V
#$ -m beas
#$ -M david@wilkox.org
#$ -pe smp 8
#
##################
perl autovalidate.pl;
perl getmean.pl < ../autovalidate_output/raw.csv > ../autovalidate_output/mean.csv
