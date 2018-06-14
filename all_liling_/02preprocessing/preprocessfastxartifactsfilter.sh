#!/bin/sh -login
#$ -cwd
#$ -V

## Join standard error and standard output into 1 file output
#$ -j y
# Time-stamp: <2014-01-16 22:24:15 liling>

# Load perl module
module load perl-5.14

usage() {
    echo "usage: preprocessfastxartifactsfilter.sh <LIB> <ADAPTER> <QUALSCORE> <MINLENGTH> <MAXLENGTH> <FXARTIFACTFILTER>"
    echo "where: <LIB> is the library name"
    echo "       <ADAPTER> is the adapter to clip"
    echo "       <QUALSCORE> is the quality score encoded in fastq file"
    echo "       <MINLENGTH> is the minimum read length"
    echo "       <MAXLENGTH> is the maximum read length"
    echo "       <FXARTIFACTFILTER> is the fastx artifact filtered file to remove reads with all but three identical bases"
    echo "       written in each row provided in a tab-separated text file."
    echo "This script does the job of running the fastx artifacts filter on the preprocessed fasta file"
    echo "from Illumina sequencing with a list of variables."
}

# Minimal argument checking
if [ $# -lt 6 ]; then
    usage
    exit
fi

# Set variables
LIB=$1        #library name
ADAPTER=$2    #adapter sequence to clip
QUALSCORE=$3  #fastq quality score according to fastQC, Q33 or Q64
MINLENGTH=$4  #18
MAXLENGTH=$5  #30
FXARTIFACTFILTER=$6  #fxaf

LENGTHFILTERFASTA=l${MINLENGTH}t${MAXLENGTH}nt

echo "Start with ${LIB}"
date

echo "print variables"
echo "<LIB>"
echo $1
echo "<ADAPTER>"
echo $2
echo "<QUALSCORE>"
echo $3
echo "<MINLENGTH>"
echo $4
echo "<MAXLENGTH>"
echo $5
echo "<FXARTIFACTFILTER>"
echo $6

echo "Do fastx artifacts filter to remove reads with all but three identical bases"
fastx_artifacts_filter -${QUALSCORE} -v -i ${LIB}clip${LENGTHFILTERFASTA}colid.fasta \
> ${LIB}clip${LENGTHFILTERFASTA}${FXARTIFACTFILTER}colid.fasta

echo "head check ${LIB}clip${LENGTHFILTERFASTA}${FXARTIFACTFILTER}colid.fasta"
head ${LIB}clip${LENGTHFILTERFASTA}${FXARTIFACTFILTER}colid.fasta

echo "count number of collapsed length filtered sequences after fastx artifacts filter"
grep -c ${LIB} ${LIB}clip${LENGTHFILTERFASTA}${FXARTIFACTFILTER}colid.fasta

echo "drawing histogram after length filter and fastx artifacts filter"
fasta_clipping_histogram2.pl ${LIB}clip${LENGTHFILTERFASTA}${FXARTIFACTFILTER}colid.fasta ${LIB}clip${LENGTHFILTERFASTA}${FXARTIFACTFILTER}colid.png

echo "Done with ${LIB}"
date
