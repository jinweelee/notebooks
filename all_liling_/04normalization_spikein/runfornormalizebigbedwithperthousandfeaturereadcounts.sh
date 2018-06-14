#!/bin/bash -login
#$ -cwd
#$ -V

## Join standard error and standard output into 1 file output
#$ -j y
# Time-stamp: <2015-12-10 15:53:59 liling>

## exit when there is an error
set -e

## set locale as C
export LC_ALL="C"

## Load perl module
module load perl-5.14

#list.txt contains <LIB> <FEATURENORM> <FEATURENORMCOUNT> <COUNTTYPE> <ASSEMBLY>

date

getlistfornormalizebigbedwithperthousandfeaturereadcounts.sh featurenorm20151210run1.txt

date
