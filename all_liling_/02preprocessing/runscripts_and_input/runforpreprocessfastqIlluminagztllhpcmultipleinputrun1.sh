#!/bin/sh -login
#$ -cwd
#$ -V

## Join standard error and standard output into 1 file output
#$ -j y
# Time-stamp: <2014-02-17 12:24:12 liling>

# Load perl module
module load perl-5.14

#list.txt contains <LIB> <ADAPTER> <QUALSCORE> <MINLENGTH> <MAXLENGTH> <FXARTIFACTFILTER>

date

getlistforpreprocessfastqgztllhpc.sh preprocess20151203run1.txt

date
