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
###################
perl thesis.pl 2>&1 > autovalidate.log
