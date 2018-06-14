#!/bin/sh -login
#$ -cwd
#$ -V
#$ -pe smp 2

## Join standard error and standard output into 1 file output
#$ -j y
# Time-stamp: <2015-05-16 18:33:35 liling>

## exit when there is an error
set -e

## set locale as C
export LC_ALL="C"

## Load perl module
module load perl-5.14

#list.txt contains <LIB> <MINLENGTH> <MAXLENGTH> <MAPINDEX> <INDEXBASE> <NTMISMATCH> <ALIGNREPORT> <BWBBTYPE> <ASSEMBLY> <FXARTIFACTFILTER>

date

getlistformapcolfastabwbbnormRPM.sh mapgenome20151203run2.txt

date
