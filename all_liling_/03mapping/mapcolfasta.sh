#!/bin/bash -login
#$ -cwd
#$ -V
#$ -pe smp 2
#$ -j y
# Time-stamp: <2015-11-09 20:07:07 liling>

## exit when there is an error
set -e

## set locale as C
export LC_ALL="C"

## load perl module
module load perl-5.14

usage () {
    echo "usage: mapcolfasta.sh <LIB> <MINLENGTH> <MAXLENGTH> <MAPINDEX> <INDEXBASE> <NTMISMATCH> <ALIGNREPORT> <FXARTIFACTFILTER>"
    echo "where: <LIB> is the library name"
    echo "       <MINLENGTH> is the minimum read length"
    echo "       <MAXLENGTH> is the maximum read length"
    echo "       <MAPINDEX> is the bowtie1 index to map"
    echo "       <INDEXBASE> is the base name for the bowtie1 index"
    echo "       <NTMISMATCH> is the number of nucleotide mismatches allowed"
    echo "       <ALIGNREPORT> is the number of alignments to report"
    echo "       <FXARTIFACTFILTER> is the fastx artifact filtered file to remove reads with all but three identical bases"
    echo "       written in each row provided in a tab-separated text file."
    echo "This script does the job of mapping collapsed fasta files"
    echo "to bowtie1 index with a list of variables,"
    echo "fixing bowtie output format, making bed files and"
    echo "calculating normalized mapped reads to genomic hits, RPM values."
    echo "Location to fasta file and bowtie indexes hardcoded in script."
}

# Minimal argument checking
if [ $# -lt 7 ]; then
    usage
    exit
fi

# Set variables
LIB=$1
MINLENGTH=$2  #18
MAXLENGTH=$3  #30
MAPINDEX=$4   #mm10onlychr
INDEXBASE=$5  #mm10onlychr
NTMISMATCH=$6 #0
ALIGNREPORT=$7 #a
FXARTIFACTFILTER=$8 #fxaf

LENGTHFILTERFASTA=l${MINLENGTH}t${MAXLENGTH}nt

BOWTIE_INDEXES=/data/OkamuraLab/local/bowtieindexes/

INPUTMAPFASTAFILE=/data/OkamuraLab/processedlibraries/clip${LENGTHFILTERFASTA}${FXARTIFACTFILTER}colidfasta/${LIB}clip${LENGTHFILTERFASTA}${FXARTIFACTFILTER}colid.fasta

#NSLOTS=1   # see pe smp above

echo "Start with ${LIB}"
date

echo "print variables"
echo "<LIB>"
echo $1
echo "<MINLENGTH>"
echo $2
echo "<MAXLENGTH>"
echo $3
echo "<MAPINDEX>"
echo $4
echo "<INDEXBASE>"
echo $5
echo "<NTMISMATCH>"
echo $6
echo "<ALIGNREPORT>"
echo $7
echo "<FXARTIFACTFILTER>"
echo $8

echo "check locale setting"
locale

echo "This is running mapcolfasta script from $PWD"

printf "\n"

echo "head check ${INPUTMAPFASTAFILE} file"
head ${INPUTMAPFASTAFILE}

printf "\n"

echo "running bowtie, ${MAPINDEX} index, fasta input, \
report ${NTMISMATCH} mismatches, ${ALIGNREPORT} alignments, bowtie output format"
bowtie -f --un ${LIB}unmapcol${INDEXBASE}.fasta -${ALIGNREPORT} -v ${NTMISMATCH} -p ${NSLOTS} -t ${MAPINDEX} ${INPUTMAPFASTAFILE} \
> ${LIB}mapcol${INDEXBASE}.bowtie.txt

echo "parameters used"
echo "bowtie -f --un ${LIB}unmapcol${INDEXBASE}.fasta -${ALIGNREPORT} -v ${NTMISMATCH} -p ${NSLOTS} -t ${MAPINDEX}"

echo "done with mapping"

echo "head check ${LIB}mapcol${INDEXBASE}.bowtie.txt"
head ${LIB}mapcol${INDEXBASE}.bowtie.txt

echo "count how many reads were mapped for ${LIB}mapcol${INDEXBASE}.bowtie.txt"
MAP=$(awk '{OFS=FS="\t"} {a[$1]++} END {for (b in a) {print b}}' ${LIB}mapcol${INDEXBASE}.bowtie.txt | \
awk 'BEGIN {OFS=FS="\t"} {split($1,idm,"_"); print idm[1],idm[2],idm[3]}' | awk '{SUM+=$3} END {if (SUM > 0) print SUM; else print "0"}')

echo $MAP

echo "count how many unmapped reads in ${LIB}unmapcol${INDEXBASE}.fasta"
UNMAP=$(fasta_formatter -i ${LIB}unmapcol${INDEXBASE}.fasta -t | \
awk 'BEGIN {OFS=FS="\t"} {split($1,idu,"_"); print idu[1],idu[2],idu[3]}' | awk '{SUM+=$3} END {if (SUM > 0) print SUM; else print "0"}')

echo $UNMAP

echo "sum of map and unmap reads"
SUMREAD=$(($MAP+$UNMAP))
echo $SUMREAD

echo "percentage of mapped reads with ${NTMISMATCH} mismatches"
PER=$(awk 'BEGIN {printf ("%0.2f", '$MAP' / '$SUMREAD' * 100)}')
echo $PER

printf "\n"

echo "append fastaCollapseID to strand orientation"
echo "add 1 to the alignment column in bowtie standard output format"
echo "do regular sort on first column"

awk 'BEGIN {OFS=FS="\t"} {print $1,$2,$3,$4,$5,$6,$7+1}' ${LIB}mapcol${INDEXBASE}.bowtie.txt > ${LIB}temp01.txt
sort -t $'\t' -k1,1 ${LIB}temp01.txt > ${LIB}temp1.txt

echo "list uniques for fastaCollapseID and strand by removing duplicates, leaving top row"
awk '!x[$1,$2]++' ${LIB}temp1.txt > ${LIB}temp2.txt

echo "sum alignmentColumn+1 if have same fastaCollapseID"
awk '{OFS=FS="\t"} {a[$1]+=$7} END {for (i in a) print i,a[i]}' ${LIB}temp2.txt > ${LIB}temp02.txt
sort -t $'\t' -k1,1 ${LIB}temp02.txt > ${LIB}temp3.txt  #This is the correct genomic hits count.

echo "convert mapping fasta file to tab format"
fasta_formatter -i ${INPUTMAPFASTAFILE} -t > ${LIB}temp03.txt
sort -t $'\t' -k1,1 ${LIB}temp03.txt > ${LIB}clip${LENGTHFILTERFASTA}${FXARTIFACTFILTER}colidtab.txt

echo "do join to match mapped id with mapping fasta file"
join -t $'\t' ${LIB}temp3.txt ${LIB}clip${LENGTHFILTERFASTA}${FXARTIFACTFILTER}colidtab.txt > ${LIB}temp4.txt

echo "do join to match bowtie standard output alignment+1"
join -t $'\t' ${LIB}temp4.txt ${LIB}temp1.txt > ${LIB}mapcol${INDEXBASE}.bowtie.fix.txt

echo "head check ${LIB}mapcol${INDEXBASE}.bowtie.fix.txt"
head ${LIB}mapcol${INDEXBASE}.bowtie.fix.txt

printf "\n"

echo "adjust to near bed format"
awk 'BEGIN {OFS=FS="\t"} {split($1,id,"_"); print $5,$6,$6+length($7),id[1],$3,id[3],$2,id[3],$4}' ${LIB}mapcol${INDEXBASE}.bowtie.fix.txt > ${LIB}temp04.txt
sort -t $'\t' -k1,1 -k2,2n ${LIB}temp04.txt > ${LIB}mapcol${INDEXBASE}.bed

echo "head check ${LIB}mapcol${INDEXBASE}.bed"
echo "chr start stop lib sequence readCounts genomicHits readCounts strand"
head ${LIB}mapcol${INDEXBASE}.bed

echo "make bed file with ReadCounts for display"
awk 'BEGIN {OFS=FS="\t"} {print $1,$2,$3,":"$4":"$5":"$6":"$7":"(length($5))"nt",$8,$9}' ${LIB}mapcol${INDEXBASE}.bed | \
sort-bed --max-mem 4G --tmpdir $PWD - > ${LIB}mapcol${INDEXBASE}ID.bed

echo "head check ${LIB}mapcol${INDEXBASE}ID.bed"
echo "<CHR> <START> <STOP> <:${LIB}:sequence:countvalue:mapHit:lengthnt> <COUNTVALUE> <STRAND>"
head ${LIB}mapcol${INDEXBASE}ID.bed

printf "\n"

echo "calculate normalized read counts by dividing with genomic hits count"
echo "round up to 4 decimal places"
awk 'BEGIN {OFS=FS="\t"} {print $1,$2,$3,$4,$5,sprintf("%.4f",($6/$7)),$7,sprintf("%.4f",($6/$7)),$9}' ${LIB}mapcol${INDEXBASE}.bed \
> ${LIB}mapcol${INDEXBASE}normReads.bed

echo "head check ${LIB}mapcol${INDEXBASE}normReads.bed"
head ${LIB}mapcol${INDEXBASE}normReads.bed

echo "make bed file with NormReads for display"
awk 'BEGIN {OFS=FS="\t"} {print $1,$2,$3,":"$4":"$5":"$6":"$7":"(length($5))"nt",$8,$9}' ${LIB}mapcol${INDEXBASE}normReads.bed | \
sort-bed --max-mem 4G --tmpdir $PWD - > ${LIB}mapcol${INDEXBASE}normReadID.bed

echo "head check ${LIB}mapcol${INDEXBASE}normReadID.bed"
echo "<CHR> <START> <STOP> <:${LIB}:sequence:countvalue:mapHit:lengthnt> <COUNTVALUE> <STRAND>"
head ${LIB}mapcol${INDEXBASE}normReadID.bed

printf "\n"

echo "calculate normalized reads per million from normalized reads divided by mapped library size times 10^6"
echo "rounding up normalized reads per million to four decimal places"
awk 'BEGIN {OFS=FS="\t"} {print $1,$2,$3,$4,$5,sprintf("%.4f",(($6/'$MAP')*1000000)),$7,sprintf("%.4f",(($8/'$MAP')*1000000)),$9}' ${LIB}mapcol${INDEXBASE}normReads.bed \
> ${LIB}mapcol${INDEXBASE}normRPMdeci.bed

echo "head check ${LIB}mapcol${INDEXBASE}normRPMdeci.bed"
head ${LIB}mapcol${INDEXBASE}normRPMdeci.bed

echo "reformat identifiers to :${LIB}:sequence:NormReadsRPM:mapHit:lengthnt"
awk 'BEGIN {OFS=FS="\t"} {print $1,$2,$3,":"$4":"$5":"$6":"$7":"(length($5))"nt",$8,$9}' ${LIB}mapcol${INDEXBASE}normRPMdeci.bed | \
sort-bed --max-mem 4G --tmpdir $PWD - > ${LIB}mapcol${INDEXBASE}normRPMdeciID.bed

echo "head check ${LIB}mapcol${INDEXBASE}normRPMdeciID.bed"
echo "<CHR> <START> <STOP> <:${LIB}:sequence:countvalue:mapHit:lengthnt> <COUNTVALUE> <STRAND>"
head ${LIB}mapcol${INDEXBASE}normRPMdeciID.bed

printf "\n"

echo "remove unneeded temp files"
rm ${LIB}temp01.txt
rm ${LIB}temp02.txt
rm ${LIB}temp03.txt
rm ${LIB}temp04.txt
rm ${LIB}temp1.txt
rm ${LIB}temp2.txt
rm ${LIB}temp3.txt
rm ${LIB}temp4.txt
rm ${LIB}clip${LENGTHFILTERFASTA}${FXARTIFACTFILTER}colidtab.txt
rm ${LIB}mapcol${INDEXBASE}normRPMdeci.bed
rm ${LIB}unmapcol${INDEXBASE}.fasta
rm ${LIB}mapcol${INDEXBASE}normReads.bed
rm ${LIB}mapcol${INDEXBASE}.bed

## comment out when necessary
echo "remove unnecessary files"
echo "remove ${LIB}mapcol${INDEXBASE}.bowtie.txt"
rm ${LIB}mapcol${INDEXBASE}.bowtie.txt

echo "remove ${LIB}mapcol${INDEXBASE}.bowtie.fix.txt"
rm ${LIB}mapcol${INDEXBASE}.bowtie.fix.txt

printf "\n"

echo "This done running mapcolfasta script from $PWD"

echo "Done with ${LIB}"
date
