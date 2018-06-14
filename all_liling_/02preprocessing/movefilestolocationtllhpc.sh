#!/bin/sh -login
#$ -cwd
#$ -V

## Join standard error and standard output into 1 file output
#$ -j y
# Time-stamp: <2014-01-16 22:22:35 liling>

# Load perl module
module load perl-5.14

usage() {
    echo "usage: movefilestolocationtllhpc.sh <LIB> <MINLENGTH> <MAXLENGTH> <FXARTIFACTFILTER>"
    echo "where: <LIB> is the library name"
    echo "       <MINLENGTH> is the minimum read length"
    echo "       <MAXLENGTH> is the maximum read length"
    echo "       <FXARTIFACTFILTER> is the fastx artifact filtered file to remove reads with all but three identical bases"
    echo "       written in each row provided in a tab-separated text file."
    echo "       Run this script in the same directory"
    echo "       as the preprocessing files."
    echo "This script does the job of moving files after preprocessing"
    echo "into the respective locations prior to mapping"
    echo "and making alignment files."
}

# Minimal argument checking

if [ $# -lt 3 ]; then
    usage
    exit
fi

# Set variables
LIB=$1        #library name
MINLENGTH=$2  #18
MAXLENGTH=$3  #30
FXARTIFACTFILTER=$4   #fxaf

LENGTHFILTERFASTA=l${MINLENGTH}t${MAXLENGTH}nt


echo "Start with ${LIB}"
date

echo "print variables"
echo "<LIB>"
echo $1
echo "<MINLENGTH>"
echo $2
echo "<MAXLENGTH>"
echo $3
echo "<FXARTIFACTFILTER>"
echo $4


echo "run this script within the same directory as the preprocessing files"
echo "check if the destination directory is present"
echo "if not, make destination directory"

if [ -d /data/OkamuraLab/processedlibraries/clip${LENGTHFILTERFASTA}${FXARTIFACTFILTER}colidfasta/ ]
then 
    echo "Directory present /data/OkamuraLab/processedlibraries/clip${LENGTHFILTERFASTA}${FXARTIFACTFILTER}colidfasta/"
else
    echo "Making /data/OkamuraLab/processedlibraries/clip${LENGTHFILTERFASTA}${FXARTIFACTFILTER}colidfasta/"
    mkdir -p /data/OkamuraLab/processedlibraries/clip${LENGTHFILTERFASTA}${FXARTIFACTFILTER}colidfasta/
fi

if [ -d /data/OkamuraLab/processedlibraries/clip${LENGTHFILTERFASTA}colidfasta/ ]
then 
    echo "Directory present /data/OkamuraLab/processedlibraries/clip${LENGTHFILTERFASTA}colidfasta/"
else
    echo "Making /data/OkamuraLab/processedlibraries/clip${LENGTHFILTERFASTA}colidfasta/"
    mkdir -p /data/OkamuraLab/processedlibraries/clip${LENGTHFILTERFASTA}colidfasta/
fi

if [ -d /data/OkamuraLab/processedlibraries/clipcolidfasta/ ]
then 
    echo "Directory present /data/OkamuraLab/processedlibraries/clipcolidfasta/"
else
    echo "Making /data/OkamuraLab/processedlibraries/clipcolidfasta/"
    mkdir -p /data/OkamuraLab/processedlibraries/clipcolidfasta/
fi

echo "move fasta files after preprocessing to the correct locations"
echo "move ${LIB}clip${LENGTHFILTERFASTA}${FXARTIFACTFILTER}colid.fasta"
mv ${LIB}clip${LENGTHFILTERFASTA}${FXARTIFACTFILTER}colid.fasta /data/OkamuraLab/processedlibraries/clip${LENGTHFILTERFASTA}${FXARTIFACTFILTER}colidfasta/

echo "move ${LIB}clip${LENGTHFILTERFASTA}colid.fasta"
mv ${LIB}clip${LENGTHFILTERFASTA}colid.fasta /data/OkamuraLab/processedlibraries/clip${LENGTHFILTERFASTA}colidfasta/

echo "move ${LIB}clipcolid.fasta"
mv ${LIB}clipcolid.fasta /data/OkamuraLab/processedlibraries/clipcolidfasta/

echo "remove ${LIB}clip.fasta"
rm ${LIB}clip.fasta

echo "remove ${LIB}.fastq"
rm ${LIB}.fastq

echo "remove ${LIB}clipcolid.fasta.tab"
rm ${LIB}clipcolid.fasta.tab

echo "Done with ${LIB}"
date
