#!/bin/bash -login
#$ -cwd
#$ -V

## Join standard error and standard output into 1 file output
#$ -j y
# Time-stamp: <2015-12-10 15:54:10 liling>

## exit when there is an error
set -e

## set locale as C
export LC_ALL="C"

## load perl module
module load perl-5.14

usage() {
    echo "usage: getlistfornormalizebigbedwithperthousandfeaturereadcounts.sh <DATAPROCESSLIST>"
    echo "where: <DATAPROCESSLIST> is a text file with each row for a library"
    echo "       to process with <LIB> <FEATURENORM> <FEATURENORMCOUNT> <COUNTTYPE> <ASSEMBLY>"
    echo "       in which"
    echo "       <LIB> is the library name"
    echo "       <FEATURENORM> is the feature to be normalized against"
    echo "       <FEATURENORMCOUNT> is the count value for the feature"
    echo "       <COUNTTYPE> is the type of count value file"
    echo "       <ASSEMBLY> is the genome assembly"
    echo "       written in each row provided in a tab-separated text file."
    echo "This script does the job of fetching the list of variables to pass to"
    echo "the script normalizebigbedwithperthousandfeaturereadcounts."
}


# Minimal argument checking

if [ $# -lt 1 ]; then
    usage
    exit
fi

# Set variable for input file
DATAPROCESSLIST=$1

echo "Start"
date

while read -r LINE
do
    LIB=$(echo "$LINE" | awk 'BEGIN {OFS=FS="\t"} {print $1}')
    FEATURENORM=$(echo "$LINE" | awk 'BEGIN {OFS=FS="\t"} {print $2}')
    FEATURENORMCOUNT=$(echo "$LINE" | awk 'BEGIN {OFS=FS="\t"} {print $3}')
    COUNTTYPE=$(echo "$LINE" | awk 'BEGIN {OFS=FS="\t"} {print $4}')
    ASSEMBLY=$(echo "$LINE" | awk 'BEGIN {OFS=FS="\t"} {print $5}')

    echo "variables to process"
    echo $LIB $FEATURENORM $FEATURENORMCOUNT $COUNTTYPE $ASSEMBLY

    printf "\n"

    echo "run normalizebigbedwithperthousandfeaturereadcounts"
    normalizebigbedwithperthousandfeaturereadcounts.sh $LIB $FEATURENORM $FEATURENORMCOUNT $COUNTTYPE $ASSEMBLY

    printf "\n"

done < $1

echo "Done"
date
