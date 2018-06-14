#!/bin/sh -login
#$ -cwd
#$ -V

## Join standard error and standard output into 1 file output
#$ -j y
# Time-stamp: <2014-02-17 12:20:32 liling>

# Load perl module
module load perl-5.14

usage() {
    echo "usage: preprocessfastqIlluminatllhpcgzmultipleinput.sh <LIB> <ADAPTER> <QUALSCORE> <MINLENGTH> <MAXLENGTH>"
    echo "where: <LIB> is the library name"
    echo "       <ADAPTER> is the adapter to clip"
    echo "       <QUALSCORE> is the quality score encoded in fastq file"
    echo "       <MINLENGTH> is the minimum read length"
    echo "       <MAXLENGTH> is the maximum read length"
    echo "       written in each row provided in a tab-separated text file."
    echo "This script does the job of preprocessing fastq files"
    echo "from Illumina sequencing with a list of variables."
    echo "Fastq file is gunzip, clipped, converted to fasta, length filtered,"
    echo "and collapsed into unique reads with count embeded in"
    echo "fasta identifier."
}

# Minimal argument checking
if [ $# -lt 5 ]; then
    usage
    exit
fi

# Set variables
LIB=$1        #library name
ADAPTER=$2    #adapter sequence to clip
QUALSCORE=$3  #fastq quality score according to fastQC, Q33 or Q64
MINLENGTH=$4  #18
MAXLENGTH=$5  #30

LENGTHFILTERFASTA=l${MINLENGTH}t${MAXLENGTH}nt


###############################################################################################
#ADAPTER=TGGAATTCTCGGGTGCCAAGG   # This is for the new Hiseq run 20130917 data collect.
#                                # and for 20131213 data collect.
## CTGTAGGCACCATCAATC     #Lai modencode libraries
## ATCTCGTATGCCGTCTTCTGCTTG    #3prime adaptor to clip from Illumina v1.5 Small RNA kits
## TCGTATGCCGTCTTCTGCTTG    #3prime adaptor to clip from Illumina v1 Small RNA kits
################################################################################################

# Script for processing steps from qseq to fastq format ready for mapping.
# Include option for DDBJ fastq file.
# Uncomment options as needed

# Note that changed name for fasta_clipping_histogram.pl to fasta_clipping_histogram2.pl
# to differentiate between the corrected version with the perl env.

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


echo "uncompress ${LIB}fastq.gz file"
gunzip -c ${LIB}.fastq.gz > ${LIB}.fastq

echo "head check file format ${LIB}"
#head ${LIB}.txt     # for qseq txt 

head ${LIB}.fastq    # for SRR fastq or DDBJ download, double-check Illumina version


###### Convert qseq format to fastq format. qseq 11th field signal purity filter, 
#  pass 1, fail 0.
#echo "converting qseq to fastq format, keeping reads that do not pass quality \
# control filter at 11th field"
#qseq2fastq.pl ${LIB}.txt > ${LIB}.fastq

# Head check file conversion
#echo "head check file conversion ${LIB}"
#head ${LIB}.fastq
#
######


echo "Processing files to clip, convert fasta, identifiers."
echo "clipping fasta file, discard sequences with no adaptor,"
echo "discard sequences with 'N's, discard sequences shorter than 5 nt..."
echo "converting fastq to fasta, rename identifiers to numbers..."
echo "changing read identifiers to include library information..."

##### Check if qseq, Illumina fastq, or Sanger fastq (download from DDBJ, SRA) #####
##### Check if the adaptor sequence is still the same or changed #####
## qseq converted to Illumina fastq uses -Q64 which is default for fastx_toolkit.
## if Sanger fastq, use -Q33 for all fastx_toolkit tools.
## clip minimum length default 5 nt.

fastx_clipper -${QUALSCORE} -a ${ADAPTER} -c -v -i ${LIB}.fastq | fastq_to_fasta -v -${QUALSCORE} -r | \
sed "s/^>/>${LIB}_/" > ${LIB}clip.fasta


echo "head check"
head ${LIB}clip.fasta


############# ORIGINAL COMMANDS TO DO LENGTH FILTER FIRST #################################
# Length filter with length_cutoff.pl from "NGS tools from the novice". 
#  Sequences that match length filter output as 'sequences_ok.fas'.
#echo "Length filter for min ${MINLENGTH}nt and max ${MAXLENGTH}nt."
#length_cutoff.pl -i ${LIB}clip.fasta -min ${MINLENGTH} -max ${MAXLENGTH} -0
#
#echo "rename filtered file"
#mv sequences_ok.fas ${LIB}clip${LENGTHFILTERFASTA}.fasta
#
#echo "remove too short and too long fasta files"
#rm sequences_too_long.fas
#rm sequences_too_short.fas
###########################################################################################

echo "Do collapsing of fasta file first to unique reads with collapsed counts in identifier"

echo "collapsing clipped fasta file"
echo "rename identifiers for ${LIB}"
fastx_collapser -v -i ${LIB}clip.fasta | \
sed "s/-/_/g ; s/^>/>${LIB}_/" > ${LIB}clipcolid.fasta

echo "head check collapsed file"
head ${LIB}clipcolid.fasta


echo "count numbers of collapsed sequences"
grep -c ${LIB} ${LIB}clipcolid.fasta

echo "Convert fasta format to tabular format"
echo "use information from header for collapsed sequence read count"
fasta_formatter -i ${LIB}clipcolid.fasta -t | \
awk 'BEGIN {OFS="\t"; FS="_|\t"} {print $1,$2,$3,$4,length($4)}' > ${LIB}clipcolidtemp01.txt

echo "sort by highest count to lowest"
sort -k3,3nr ${LIB}clipcolidtemp01.txt > ${LIB}clipcolid.fasta.tab

echo "head check ${LIB}clipcolid.fasta.tab"
echo "LIB SERIALNUMBER COUNT SEQUENCE LENGTH"
head -20 ${LIB}clipcolid.fasta.tab
wc -l ${LIB}clipcolid.fasta.tab

echo "remove temp files"
rm ${LIB}clipcolidtemp01.txt


echo "Calculate length counts"
awk 'BEGIN {OFS=FS="\t"} {a[$5]+=$3} END {for (i in a) print i,a[i]}' ${LIB}clipcolid.fasta.tab \
> ${LIB}clipcolidtemp02.txt

echo "sort by smallest length to highest"
sort -k1,1n ${LIB}clipcolidtemp02.txt > ${LIB}clipcolid.fasta.lengthcounts.txt

echo "head check ${LIB}clipcolid.fasta.lengthcounts.txt"
echo "LENGTH COUNT"
head -100 ${LIB}clipcolid.fasta.lengthcounts.txt

echo "remove temp files"
rm ${LIB}clipcolidtemp02.txt


echo "Length filter for min ${MINLENGTH}nt and max ${MAXLENGTH}nt"
awk '($5 >= '$MINLENGTH') && ($5 <= '$MAXLENGTH')' ${LIB}clipcolid.fasta.tab > ${LIB}clipcolid${LENGTHFILTERFASTA}temp03.txt

head ${LIB}clipcolid${LENGTHFILTERFASTA}temp03.txt

echo "convert to fasta format"
awk 'BEGIN {OFS=FS="\t"} {print $1"_"$2"_"$3,$4}' ${LIB}clipcolid${LENGTHFILTERFASTA}temp03.txt | \
awk '{print ">"$1"\n"$2}' > ${LIB}clip${LENGTHFILTERFASTA}colid.fasta

head ${LIB}clip${LENGTHFILTERFASTA}colid.fasta

echo "count number of collapsed length filtered sequences"
grep -c ${LIB} ${LIB}clip${LENGTHFILTERFASTA}colid.fasta

echo "Count total reads in ${LIB}clipcolid.fasta.tab"
awk '{SUM+=$3} END {if (SUM > 0) print SUM; else print "0"}' ${LIB}clipcolid.fasta.tab

echo "Count reads min ${MINLENGTH}nt max ${MAXLENGTH}nt"
awk '{SUM+=$3} END {if (SUM > 0) print SUM; else print "0"}' ${LIB}clipcolid${LENGTHFILTERFASTA}temp03.txt

echo "Count reads shorter than ${MINLENGTH}nt"
awk '($5 < '$MINLENGTH')' ${LIB}clipcolid.fasta.tab | awk '{SUM+=$3} END {if (SUM > 0) print SUM; else print "0"}'

echo "Count reads longer than ${MAXLENGTH}nt"
awk '($5 > '$MAXLENGTH')' ${LIB}clipcolid.fasta.tab | awk '{SUM+=$3} END {if (SUM > 0) print SUM; else print "0"}'


echo "remove temp files"
rm ${LIB}clipcolid${LENGTHFILTERFASTA}temp03.txt


echo "Done with ${LIB}"
date
