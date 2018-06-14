#!/bin/sh -login
#$ -cwd
#$ -V
#$ -pe smp 2

## Join standard error and standard output into 1 file output
#$ -j y
# Time-stamp: <2014-03-24 16:07:21 liling>

# Load perl module
module load perl-5.14

usage() {
    echo "usage: getlistforlibstatmapping.sh <DATAPROCESSLIST>"
    echo "where: <DATAPROCESSLIST> is a text file with each row for a library"
    echo "       to process with <LIB> <MINLENGTH> <MAXLENGTH> <MAPINDEX> <INDEXBASE> <NTMISMATCH> <ALIGNREPORT> <BWBBTYPE> <ASSEMBLY> <FXARTIFACTFILTER>"
    echo "       in which"
    echo "       <LIB> is the library name"
    echo "       <MINLENGTH> is the minimum read length"
    echo "       <MAXLENGTH> is the maximum read length"
    echo "       <MAPINDEX> is the bowtie1 index to map"
    echo "       <INDEXBASE> is the base name for the bowtie1 index"
    echo "       <NTMISMATCH> is the number of nucleotide mismatches allowed"
    echo "       <ALIGNREPORT> is the number of alignments to report"
    echo "       <BWBBTYPE> is the type of bigwig and bigbed file to make"
    echo "       <ASSEMBLY> is the assembly name"
    echo "       <FXARTIFACTFILTER> is the fastx artifact filtered file to remove reads with all but three identical bases"
    echo "       written in each row provided in a tab-separated text file."
    echo "This script does the job of fetching the list of variables to pass to"
    echo "the script runlibstatmapping.sh for counting numbers."
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
    LIB=$(echo $LINE | awk '{print $1}')
    MINLENGTH=$(echo $LINE | awk '{print $2}')
    MAXLENGTH=$(echo $LINE | awk '{print $3}')
    MAPINDEX=$(echo $LINE | awk '{print $4}')
    INDEXBASE=$(echo $LINE | awk '{print $5}')
    NTMISMATCH=$(echo $LINE | awk '{print $6}')
    ALIGNREPORT=$(echo $LINE | awk '{print $7}')
    BWBBTYPE=$(echo $LINE | awk '{print $8}')
    ASSEMBLY=$(echo $LINE | awk '{print $9}')
    FXARTIFACTFILTER=$(echo $LINE | awk '{print $10}')
    echo "variables to process"
    echo $LIB $MINLENGTH $MAXLENGTH $MAPINDEX $INDEXBASE $NTMISMATCH $ALIGNREPORT $BWBBTYPE $ASSEMBLY $FXARTIFACTFILTER
    echo "run mapping bowtie1"
    runlibstatmapping.sh $LIB $MINLENGTH $MAXLENGTH $MAPINDEX $INDEXBASE $NTMISMATCH $ALIGNREPORT $FXARTIFACTFILTER
done < $1

echo "Done"
date
