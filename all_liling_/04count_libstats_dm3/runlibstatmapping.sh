#!/bin/sh -login
#$ -cwd
#$ -V
#$ -pe smp 2

## Join standard error and standard output into 1 file output
#$ -j y
# Time-stamp: <2014-03-25 12:16:22 liling>

# Load perl module
module load perl-5.14

usage () {
    echo "usage: runlibstatmapping.sh <LIB> <MINLENGTH> <MAXLENGTH> <MAPINDEX> <INDEXBASE> <NTMISMATCH> <ALIGNREPORT> <FXARTIFACTFILTER>"
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
    echo "with subsequent unmatched fasta files as input for the next mapping."
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
MAPINDEX=$4   #dme3excludeUextra
INDEXBASE=$5  #dm3excludeUext
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

# Set variables
INDEXBASE=dm3excludeUext
FBDMELALLGFF=dmel-all-r5.52.gff
FBDMELBASE=FBdmelr5.52
FBDMELBED=FBdmelr5.52.bed
FBEXPT1=${FBDMELBASE}.txt.exon.mRNA.merge  ## 1749 overlaps with mRNA
FBEXPT2=${FBDMELBASE}.txt.exon.ncRNA.merge ## 492 overlaps with mRNA
FBEXPT3=${FBDMELBASE}.txt.exon.rRNA.merge ## no overlaps
FBEXPT4=${FBDMELBASE}.txt.exon.tRNA.merge  ## 2 overlaps with mRNA
FBEXPT5=${FBDMELBASE}.txt.exon.snRNA.merge  ## no overlaps
FBEXPT6=${FBDMELBASE}.txt.exon.snoRNA.merge  ## 5 overlaps with mRNA
FBEXPT7=${FBDMELBASE}.txt.exon.pseudogene.merge  ## 16 overlaps with mRNA
FBDMELBEDTYPE1=${FBDMELBED}.PT.exon ## this separates into mRNA, ncRNA, rRNA, pseudogene, tRNA, snoRNA, snRNA
FBDMELBEDTYPE2=${FBDMELBED}.PT.intron ## this separates into mRNA, ncRNA, pseudogene, tRNA
FBDMELBEDTYPE3=${FBDMELBED}.PT.threeprimeUTR
FBDMELBEDTYPE4=${FBDMELBED}.PT.fiveprimeUTR
FBDMELBEDTYPE5=${FBDMELBED}.PT.transposableelement

ANNOTRMSK1=LTR
ANNOTRMSK2=LINE
ANNOTRMSK3=DNA
ANNOTRMSK4=Satellite
ANNOTRMSK5=Low_complexity
ANNOTRMSK6=rRNA
ANNOTRMSK7=RNA
ANNOTRMSK8=RC
ANNOTRMSK9=Simple_repeat
ANNOTRMSK10=Other
ANNOTRMSK11=Unknown
ANNOTRMSK12=ARTEFACT

FBANNOT3=${FBDMELBEDTYPE1}.rRNA.merge ## exon
FBANNOT5=${FBDMELBEDTYPE1}.tRNA.merge ## exon
FBANNOT7=${FBDMELBEDTYPE1}.snRNA.merge ## exon
FBANNOT6=${FBDMELBEDTYPE1}.snoRNA.merge ## exon
FBANNOT16=${FBDMELBEDTYPE5}.transposableelement.merge ## This is further annotation for transposable_elements.

FBOVERLAPANNOT1=overlap${FBEXPT1}.${FBEXPT1} ## 1749 overlaps with mRNA.mRNA
FBOVERLAPANNOT2=overlap${FBEXPT1}.${FBEXPT2} ## 492 overlaps with mRNA.ncRNA
FBOVERLAPANNOT4=overlap${FBEXPT1}.${FBEXPT4} ## 2 overlaps with mRNA.tRNA
FBOVERLAPANNOT6=overlap${FBEXPT1}.${FBEXPT6} ## 5 overlaps with mRNA.snoRNA
FBOVERLAPANNOT7=overlap${FBEXPT1}.${FBEXPT7} ## 16 overlaps with mRNA.pseudogene
FBANNOT2=${FBDMELBEDTYPE1}.ncRNA.merge ## exon
FBANNOT4=${FBDMELBEDTYPE1}.pseudogene.merge ## exon


## Set bowtie indexes

# STEP 1 : background RNAs
RMSKINDEX6=${INDEXBASE}${ANNOTRMSK6} #repeatmasker local run dm3excludeUextrRNA
FBINDEX3=${INDEXBASE}${FBANNOT3} # exon rRNA
DMERDNAINDEX=chr1dmerdna # dme ribosomal DNA, NCBI #index
RMSKINDEX7=${INDEXBASE}${ANNOTRMSK7}   #repeatmasker local run dm3excludeUextRNA
FBINDEX5=${INDEXBASE}${FBANNOT5} # exon tRNA
FBINDEX7=${INDEXBASE}${FBANNOT7} # exon snRNA
FBINDEX6=${INDEXBASE}${FBANNOT6} # exon snoRNA

INDEX1=${RMSKINDEX6}
INDEX2=${FBINDEX3}
INDEX3=${DMERDNAINDEX}
INDEX4=${RMSKINDEX7}
INDEX5=${FBINDEX5}
INDEX6=${FBINDEX7}
INDEX7=${FBINDEX6}


# STEP 2 : mirprecursor
MIRPREINDEX=dmestemloopmirbasev20 #index

INDEX8=${MIRPREINDEX}


# STEP 3 : known hairpin small RNAs  ## Okamura2008nature # INDEX made
HAIRPINRNAINDEX=knownhprna #index

INDEX9=${HAIRPINRNAINDEX}


# STEP 4 : transposons, transposable elements, repetitive regions
RMSKINDEX1=${INDEXBASE}${ANNOTRMSK1}  #repeatmasker local run LTR
RMSKINDEX2=${INDEXBASE}${ANNOTRMSK2}  #repeatmasker local run LINE
RMSKINDEX3=${INDEXBASE}${ANNOTRMSK3}  #repeatmasker local run DNA
RMSKINDEX4=${INDEXBASE}${ANNOTRMSK4}  #repeatmasker local run Satellite
RMSKINDEX5=${INDEXBASE}${ANNOTRMSK5}  #repeatmasker local run Low_complexity
RMSKINDEX8=${INDEXBASE}${ANNOTRMSK8}  #repeatmasker local run RC
RMSKINDEX9=${INDEXBASE}${ANNOTRMSK9}  #repeatmasker local run Simple_repeat
RMSKINDEX10=${INDEXBASE}${ANNOTRMSK10}  #repeatmasker local run Other
RMSKINDEX11=${INDEXBASE}${ANNOTRMSK11}  #repeatmasker local run Unknown
FBINDEX16=${INDEXBASE}${FBANNOT16} # transposable elements
RMSKINDEX12=${INDEXBASE}${ANNOTRMSK12}  #repeatmasker local run ARTEFACT

INDEX10=${RMSKINDEX1}
INDEX11=${RMSKINDEX2}
INDEX12=${RMSKINDEX3}
INDEX13=${RMSKINDEX4}
INDEX14=${RMSKINDEX5}
INDEX15=${RMSKINDEX8}
INDEX16=${RMSKINDEX9}
INDEX17=${RMSKINDEX10}
INDEX18=${RMSKINDEX11}
INDEX19=${FBINDEX16}
INDEX20=${RMSKINDEX12}


# STEP 5 : published small RNA including endo-siRNA loci
LOCIPUBINDEX1=${INDEXBASE}okamura2008nsmb    ## use only the 117 annotated 3' cis-NAT loci list, overlapping gene pairs, exons.
LOCIPUBINDEX2=${INDEXBASE}czech2008nature    ## AGO2-associated small RNA 20-22nt, unique map. Note, this may have non cis-NAT loci.
LOCIPUBINDEX3=${INDEXBASE}kawamura2008nature ## AGO2-associated small RNA. Note, this may have non cis-NAT loci.
LOCIPUBINDEX4=${INDEXBASE}ghildiyal2008science ## siRNAs, S2, head. Take endo siRNAs list. Note, this may have non cis-NAT loci.

INDEX21=${LOCIPUBINDEX1}
INDEX22=${LOCIPUBINDEX2}
INDEX23=${LOCIPUBINDEX3}
INDEX24=${LOCIPUBINDEX4}


# index for CR14033, removed because overlaps with the intronic region of thickveins. Known to produce large numbers of siRNAs.
CR14033INDEX=${INDEXBASE}CR14033
INDEX25=${CR14033INDEX}


# STEP 6 : exons overlap mRNA and non-mRNA, CISNAT
FBOVERLAPINDEX1=${INDEXBASE}${FBOVERLAPANNOT1} ## 1749 overlaps with mRNA.mRNA
FBOVERLAPINDEX2=${INDEXBASE}${FBOVERLAPANNOT2} ## 492 overlaps with mRNA.ncRNA
FBOVERLAPINDEX4=${INDEXBASE}${FBOVERLAPANNOT4} ## 2 overlaps with mRNA.tRNA
FBOVERLAPINDEX6=${INDEXBASE}${FBOVERLAPANNOT6} ## 5 overlaps with mRNA.snoRNA
FBOVERLAPINDEX7=${INDEXBASE}${FBOVERLAPANNOT7} ## 16 overlaps with mRNA.pseudogene

INDEX26=${FBOVERLAPINDEX1}
INDEX27=${FBOVERLAPINDEX2}
INDEX28=${FBOVERLAPINDEX4}
INDEX29=${FBOVERLAPINDEX6}
INDEX30=${FBOVERLAPINDEX7}


###### PREVIOUS for calculating numbers
# STEP 7 : exon ncRNA, pseudogene
#
#FBINDEX2=${INDEXBASE}${FBANNOT2} # exon ncRNA
#FBINDEX4=${INDEXBASE}${FBANNOT4} # exon pseudogene
#
#INDEX30=${FBINDEX2}
#INDEX31=${FBINDEX4}
######


# remap unmap category index back to the dm3excludeUext
INDEX31=categoryunknownremapdm3excludeUext


# Do
echo "Start with ${LIB}"
date

echo "head check ${INPUTMAPFASTAFILE} file"
head ${INPUTMAPFASTAFILE}


# Run bowtie mapping with ${MAPINDEX}
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
MAPGENOME=$(awk '{OFS=FS="\t"} {a[$1]++} END {for (b in a) {print b}}' ${LIB}mapcol${INDEXBASE}.bowtie.txt | \
awk 'BEGIN {FS="_"; OFS="\t"} {print $1,$2,$3}' | awk '{SUM+=$3} END {if (SUM > 0) print SUM; else print "0"}')

echo $MAPGENOME

echo "count how many unmapped reads in ${LIB}unmapcol${INDEXBASE}.fasta"
UNMAPGENOME=$(fasta_formatter -i ${LIB}unmapcol${INDEXBASE}.fasta -t | \
awk 'BEGIN {FS="_"; OFS="\t"} {print $1,$2,$3}' | awk '{SUM+=$3} END {if (SUM > 0) print SUM; else print "0"}')

echo $UNMAPGENOME

echo "sum of map and unmap reads"
SUMREADGENOME=$(($MAPGENOME+$UNMAPGENOME))
echo $SUMREADGENOME

echo "percentage of mapped reads with ${NTMISMATCH} mismatches"
PERGENOME=$(awk 'BEGIN {printf ("%0.2f", '$MAPGENOME' / '$SUMREADGENOME' * 100)}')
echo $PERGENOME


# Fix bowtie mapping file
echo "append fastaCollapseID to strand orientation"
echo "add 1 to the alignment column in bowtie standard output format"
echo "do regular sort on first column"

awk 'BEGIN {OFS=FS="\t"} {print $1,$2,$3,$4,$5,$6,$7+1}' ${LIB}mapcol${INDEXBASE}.bowtie.txt > ${LIB}temp01.txt
sort -k1,1 ${LIB}temp01.txt > ${LIB}temp1.txt

echo "list uniques for fastaCollapseID and strand by removing duplicates, leaving top row"
awk '!x[$1,$2]++' ${LIB}temp1.txt > ${LIB}temp2.txt

echo "sum alignmentColumn+1 if have same fastaCollapseID"
awk '{OFS=FS="\t"} {a[$1]+=$7} END {for (i in a) print i,a[i]}' ${LIB}temp2.txt > ${LIB}temp02.txt
sort -k1,1 ${LIB}temp02.txt > ${LIB}temp3.txt  #This is the correct genomic hits count.

echo "convert mapping fasta file to tab format"
fasta_formatter -i ${INPUTMAPFASTAFILE} -t > ${LIB}temp03.txt
sort -k1,1 ${LIB}temp03.txt > ${LIB}clip${LENGTHFILTERFASTA}${FXARTIFACTFILTER}colidtab.txt

echo "do join to match mapped id with mapping fasta file"
join ${LIB}temp3.txt ${LIB}clip${LENGTHFILTERFASTA}${FXARTIFACTFILTER}colidtab.txt > ${LIB}temp4.txt

echo "do join to match bowtie standard output alignment+1"
join ${LIB}temp4.txt ${LIB}temp1.txt > ${LIB}mapcol${INDEXBASE}.bowtie.fix.txt

echo "head check ${LIB}mapcol${INDEXBASE}.bowtie.fix.txt"
head ${LIB}mapcol${INDEXBASE}.bowtie.fix.txt


# Count unique and multiple mappers
echo "separate lists into multiple mappers and unique mappers"
awk '$2 == 1' ${LIB}mapcol${INDEXBASE}.bowtie.fix.txt > uniqmap${LIB}mapcol${INDEXBASE}.bowtie.fix.txt
awk '$2 > 1' ${LIB}mapcol${INDEXBASE}.bowtie.fix.txt > multimap${LIB}mapcol${INDEXBASE}.bowtie.fix.txt

echo "count reads that are unique mappers"
UNIMAP=$(awk '!x[$1]++' uniqmap${LIB}mapcol${INDEXBASE}.bowtie.fix.txt | \
awk 'BEGIN {OFS="\t"; FS="_| "} {print $3}' | \
awk '{SUM+=$1} END {if (SUM > 0) print SUM; else print "0"}')

echo $UNIMAP

echo "count reads that are multiple mappers"
MULMAP=$(awk '!x[$1]++' multimap${LIB}mapcol${INDEXBASE}.bowtie.fix.txt | \
awk 'BEGIN {OFS="\t"; FS="_| "} {print $3}' | \
awk '{SUM+=$1} END {if (SUM > 0) print SUM; else print "0"}')

echo $MULMAP

echo "count sum of unique and multiple mappers"
SUMREADMAPPERS=$(($UNIMAP+$MULMAP))

echo $SUMREADMAPPERS

echo "complete list"
echo $UNIMAP
echo $MULMAP
echo $SUMREADMAPPERS


# convert bowtie fix file to fasta for remapping
echo "print fastaCollapseID and sequence"
echo "list uniques for fastaCollapseID removing duplicates, leaving top row"
awk 'BEGIN {OFS="\t"; FS=" "} {print $1,$3}' ${LIB}mapcol${INDEXBASE}.bowtie.fix.txt | \
awk '!x[$1,$2]++' > ${LIB}mapcol${INDEXBASE}.tab

echo "head check ${LIB}mapcol${INDEXBASE}.tab"
head ${LIB}mapcol${INDEXBASE}.tab
wc -l ${LIB}mapcol${INDEXBASE}.tab

echo "convert to fasta"
awk '{print ">"$1"\n"$2}' ${LIB}mapcol${INDEXBASE}.tab > ${LIB}mapcol${INDEXBASE}.fasta

echo "head check ${LIB}mapcol${INDEXBASE}.fasta"
head ${LIB}mapcol${INDEXBASE}.fasta
wc -l ${LIB}mapcol${INDEXBASE}.fasta


# library statistics mapping with perfect matches, all alignments reported
echo "use reads that match perfectly to ${INDEXBASE}"
echo "head check ${LIB}mapcol${INDEXBASE}.fasta"
head ${LIB}mapcol${INDEXBASE}.fasta

echo " STEP 1 background RNAs"

echo "running bowtie, ${INDEX1} index, fasta input, \
report perfect match, all alignments, bowtie output format"
bowtie -f --un ${LIB}unmapcol${INDEX1}.fasta -a -v 0 -p ${NSLOTS} -t ${INDEX1} ${LIB}mapcol${INDEXBASE}.fasta \
> ${LIB}mapcol${INDEX1}.bowtie.txt

echo "done with mapping"

echo "head check ${LIB}mapcol${INDEX1}.bowtie.txt"
head ${LIB}mapcol${INDEX1}.bowtie.txt

echo "count how many reads were mapped for ${LIB}mapcol${INDEX1}.bowtie.txt"
MAP1=$(awk '{OFS=FS="\t"} {a[$1]++} END {for (b in a) {print b}}' ${LIB}mapcol${INDEX1}.bowtie.txt | \
awk 'BEGIN {FS="_"; OFS="\t"} {print $1,$2,$3}' | awk '{SUM+=$3} END {if (SUM > 0) print SUM; else print "0"}')

echo $MAP1


echo "running bowtie, ${INDEX2} index, fasta input from previous unmatched, \
report perfect match, all alignments, bowtie output format"
bowtie -f --un ${LIB}unmapcol${INDEX2}.fasta -a -v 0 -p ${NSLOTS} -t ${INDEX2} ${LIB}unmapcol${INDEX1}.fasta \
> ${LIB}mapcol${INDEX2}.bowtie.txt

echo "done with mapping"

echo "head check ${LIB}mapcol${INDEX2}.bowtie.txt"
head ${LIB}mapcol${INDEX2}.bowtie.txt

echo "count how many reads were mapped for ${LIB}mapcol${INDEX2}.bowtie.txt"
MAP2=$(awk '{OFS=FS="\t"} {a[$1]++} END {for (b in a) {print b}}' ${LIB}mapcol${INDEX2}.bowtie.txt | \
awk 'BEGIN {FS="_"; OFS="\t"} {print $1,$2,$3}' | awk '{SUM+=$3} END {if (SUM > 0) print SUM; else print "0"}')

echo $MAP2


echo "running bowtie, ${INDEX3} index, fasta input from previous unmatched, \
report perfect match, all alignments, bowtie output format"
bowtie -f --un ${LIB}unmapcol${INDEX3}.fasta -a -v 0 -p ${NSLOTS} -t ${INDEX3} ${LIB}unmapcol${INDEX2}.fasta \
> ${LIB}mapcol${INDEX3}.bowtie.txt

echo "done with mapping"

echo "head check ${LIB}mapcol${INDEX3}.bowtie.txt"
head ${LIB}mapcol${INDEX3}.bowtie.txt

echo "count how many reads were mapped for ${LIB}mapcol${INDEX3}.bowtie.txt"
MAP3=$(awk '{OFS=FS="\t"} {a[$1]++} END {for (b in a) {print b}}' ${LIB}mapcol${INDEX3}.bowtie.txt | \
awk 'BEGIN {FS="_"; OFS="\t"} {print $1,$2,$3}' | awk '{SUM+=$3} END {if (SUM > 0) print SUM; else print "0"}')

echo $MAP3


echo "running bowtie, ${INDEX4} index, fasta input from previous unmatched, \
report perfect match, all alignments, bowtie output format"
bowtie -f --un ${LIB}unmapcol${INDEX4}.fasta -a -v 0 -p ${NSLOTS} -t ${INDEX4} ${LIB}unmapcol${INDEX3}.fasta \
> ${LIB}mapcol${INDEX4}.bowtie.txt

echo "done with mapping"

echo "head check ${LIB}mapcol${INDEX4}.bowtie.txt"
head ${LIB}mapcol${INDEX4}.bowtie.txt

echo "count how many reads were mapped for ${LIB}mapcol${INDEX4}.bowtie.txt"
MAP4=$(awk '{OFS=FS="\t"} {a[$1]++} END {for (b in a) {print b}}' ${LIB}mapcol${INDEX4}.bowtie.txt | \
awk 'BEGIN {FS="_"; OFS="\t"} {print $1,$2,$3}' | awk '{SUM+=$3} END {if (SUM > 0) print SUM; else print "0"}')

echo $MAP4


echo "running bowtie, ${INDEX5} index, fasta input from previous unmatched, \
report perfect match, all alignments, bowtie output format"
bowtie -f --un ${LIB}unmapcol${INDEX5}.fasta -a -v 0 -p ${NSLOTS} -t ${INDEX5} ${LIB}unmapcol${INDEX4}.fasta \
> ${LIB}mapcol${INDEX5}.bowtie.txt

echo "done with mapping"

echo "head check ${LIB}mapcol${INDEX5}.bowtie.txt"
head ${LIB}mapcol${INDEX5}.bowtie.txt

echo "count how many reads were mapped for ${LIB}mapcol${INDEX5}.bowtie.txt"
MAP5=$(awk '{OFS=FS="\t"} {a[$1]++} END {for (b in a) {print b}}' ${LIB}mapcol${INDEX5}.bowtie.txt | \
awk 'BEGIN {FS="_"; OFS="\t"} {print $1,$2,$3}' | awk '{SUM+=$3} END {if (SUM > 0) print SUM; else print "0"}')

echo $MAP5


echo "running bowtie, ${INDEX6} index, fasta input from previous unmatched, \
report perfect match, all alignments, bowtie output format"
bowtie -f --un ${LIB}unmapcol${INDEX6}.fasta -a -v 0 -p ${NSLOTS} -t ${INDEX6} ${LIB}unmapcol${INDEX5}.fasta \
> ${LIB}mapcol${INDEX6}.bowtie.txt

echo "done with mapping"

echo "head check ${LIB}mapcol${INDEX6}.bowtie.txt"
head ${LIB}mapcol${INDEX6}.bowtie.txt

echo "count how many reads were mapped for ${LIB}mapcol${INDEX6}.bowtie.txt"
MAP6=$(awk '{OFS=FS="\t"} {a[$1]++} END {for (b in a) {print b}}' ${LIB}mapcol${INDEX6}.bowtie.txt | \
awk 'BEGIN {FS="_"; OFS="\t"} {print $1,$2,$3}' | awk '{SUM+=$3} END {if (SUM > 0) print SUM; else print "0"}')

echo $MAP6


echo "running bowtie, ${INDEX7} index, fasta input from previous unmatched, \
report perfect match, all alignments, bowtie output format"
bowtie -f --un ${LIB}unmapcol${INDEX7}.fasta -a -v 0 -p ${NSLOTS} -t ${INDEX7} ${LIB}unmapcol${INDEX6}.fasta \
> ${LIB}mapcol${INDEX7}.bowtie.txt

echo "done with mapping"

echo "head check ${LIB}mapcol${INDEX7}.bowtie.txt"
head ${LIB}mapcol${INDEX7}.bowtie.txt

echo "count how many reads were mapped for ${LIB}mapcol${INDEX7}.bowtie.txt"
MAP7=$(awk '{OFS=FS="\t"} {a[$1]++} END {for (b in a) {print b}}' ${LIB}mapcol${INDEX7}.bowtie.txt | \
awk 'BEGIN {FS="_"; OFS="\t"} {print $1,$2,$3}' | awk '{SUM+=$3} END {if (SUM > 0) print SUM; else print "0"}')

echo $MAP7


echo " STEP 2 mirprecursor"

echo "running bowtie, ${INDEX8} index, fasta input from previous unmatched, \
report perfect match, all alignments, bowtie output format"
bowtie -f --un ${LIB}unmapcol${INDEX8}.fasta -a -v 0 -p ${NSLOTS} -t ${INDEX8} ${LIB}unmapcol${INDEX7}.fasta \
> ${LIB}mapcol${INDEX8}.bowtie.txt

echo "done with mapping"

echo "head check ${LIB}mapcol${INDEX8}.bowtie.txt"
head ${LIB}mapcol${INDEX8}.bowtie.txt

echo "count how many reads were mapped for ${LIB}mapcol${INDEX8}.bowtie.txt"
MAP8=$(awk '{OFS=FS="\t"} {a[$1]++} END {for (b in a) {print b}}' ${LIB}mapcol${INDEX8}.bowtie.txt | \
awk 'BEGIN {FS="_"; OFS="\t"} {print $1,$2,$3}' | awk '{SUM+=$3} END {if (SUM > 0) print SUM; else print "0"}')

echo $MAP8


echo " STEP 3 known hairpin small RNAs"

echo "running bowtie, ${INDEX9} index, fasta input from previous unmatched, \
report perfect match, all alignments, bowtie output format"
bowtie -f --un ${LIB}unmapcol${INDEX9}.fasta -a -v 0 -p ${NSLOTS} -t ${INDEX9} ${LIB}unmapcol${INDEX8}.fasta \
> ${LIB}mapcol${INDEX9}.bowtie.txt

echo "done with mapping"

echo "head check ${LIB}mapcol${INDEX9}.bowtie.txt"
head ${LIB}mapcol${INDEX9}.bowtie.txt

echo "count how many reads were mapped for ${LIB}mapcol${INDEX9}.bowtie.txt"
MAP9=$(awk '{OFS=FS="\t"} {a[$1]++} END {for (b in a) {print b}}' ${LIB}mapcol${INDEX9}.bowtie.txt | \
awk 'BEGIN {FS="_"; OFS="\t"} {print $1,$2,$3}' | awk '{SUM+=$3} END {if (SUM > 0) print SUM; else print "0"}')

echo $MAP9


echo " STEP 4 transposons, transposable elements, repetitive regions"

echo "running bowtie, ${INDEX10} index, fasta input from previous unmatched, \
report perfect match, all alignments, bowtie output format"
bowtie -f --un ${LIB}unmapcol${INDEX10}.fasta -a -v 0 -p ${NSLOTS} -t ${INDEX10} ${LIB}unmapcol${INDEX9}.fasta \
> ${LIB}mapcol${INDEX10}.bowtie.txt

echo "done with mapping"

echo "head check ${LIB}mapcol${INDEX10}.bowtie.txt"
head ${LIB}mapcol${INDEX10}.bowtie.txt

echo "count how many reads were mapped for ${LIB}mapcol${INDEX10}.bowtie.txt"
MAP10=$(awk '{OFS=FS="\t"} {a[$1]++} END {for (b in a) {print b}}' ${LIB}mapcol${INDEX10}.bowtie.txt | \
awk 'BEGIN {FS="_"; OFS="\t"} {print $1,$2,$3}' | awk '{SUM+=$3} END {if (SUM > 0) print SUM; else print "0"}')

echo $MAP10


echo "running bowtie, ${INDEX11} index, fasta input from previous unmatched, \
report perfect match, all alignments, bowtie output format"
bowtie -f --un ${LIB}unmapcol${INDEX11}.fasta -a -v 0 -p ${NSLOTS} -t ${INDEX11} ${LIB}unmapcol${INDEX10}.fasta \
> ${LIB}mapcol${INDEX11}.bowtie.txt

echo "done with mapping"

echo "head check ${LIB}mapcol${INDEX11}.bowtie.txt"
head ${LIB}mapcol${INDEX11}.bowtie.txt

echo "count how many reads were mapped for ${LIB}mapcol${INDEX11}.bowtie.txt"
MAP11=$(awk '{OFS=FS="\t"} {a[$1]++} END {for (b in a) {print b}}' ${LIB}mapcol${INDEX11}.bowtie.txt | \
awk 'BEGIN {FS="_"; OFS="\t"} {print $1,$2,$3}' | awk '{SUM+=$3} END {if (SUM > 0) print SUM; else print "0"}')

echo $MAP11


echo "running bowtie, ${INDEX12} index, fasta input from previous unmatched, \
report perfect match, all alignments, bowtie output format"
bowtie -f --un ${LIB}unmapcol${INDEX12}.fasta -a -v 0 -p ${NSLOTS} -t ${INDEX12} ${LIB}unmapcol${INDEX11}.fasta \
> ${LIB}mapcol${INDEX12}.bowtie.txt

echo "done with mapping"

echo "head check ${LIB}mapcol${INDEX12}.bowtie.txt"
head ${LIB}mapcol${INDEX12}.bowtie.txt

echo "count how many reads were mapped for ${LIB}mapcol${INDEX12}.bowtie.txt"
MAP12=$(awk '{OFS=FS="\t"} {a[$1]++} END {for (b in a) {print b}}' ${LIB}mapcol${INDEX12}.bowtie.txt | \
awk 'BEGIN {FS="_"; OFS="\t"} {print $1,$2,$3}' | awk '{SUM+=$3} END {if (SUM > 0) print SUM; else print "0"}')

echo $MAP12



echo "running bowtie, ${INDEX13} index, fasta input from previous unmatched, \
report perfect match, all alignments, bowtie output format"
bowtie -f --un ${LIB}unmapcol${INDEX13}.fasta -a -v 0 -p ${NSLOTS} -t ${INDEX13} ${LIB}unmapcol${INDEX12}.fasta \
> ${LIB}mapcol${INDEX13}.bowtie.txt

echo "done with mapping"

echo "head check ${LIB}mapcol${INDEX13}.bowtie.txt"
head ${LIB}mapcol${INDEX13}.bowtie.txt

echo "count how many reads were mapped for ${LIB}mapcol${INDEX13}.bowtie.txt"
MAP13=$(awk '{OFS=FS="\t"} {a[$1]++} END {for (b in a) {print b}}' ${LIB}mapcol${INDEX13}.bowtie.txt | \
awk 'BEGIN {FS="_"; OFS="\t"} {print $1,$2,$3}' | awk '{SUM+=$3} END {if (SUM > 0) print SUM; else print "0"}')

echo $MAP13



echo "running bowtie, ${INDEX14} index, fasta input from previous unmatched, \
report perfect match, all alignments, bowtie output format"
bowtie -f --un ${LIB}unmapcol${INDEX14}.fasta -a -v 0 -p ${NSLOTS} -t ${INDEX14} ${LIB}unmapcol${INDEX13}.fasta \
> ${LIB}mapcol${INDEX14}.bowtie.txt

echo "done with mapping"

echo "head check ${LIB}mapcol${INDEX14}.bowtie.txt"
head ${LIB}mapcol${INDEX14}.bowtie.txt

echo "count how many reads were mapped for ${LIB}mapcol${INDEX14}.bowtie.txt"
MAP14=$(awk '{OFS=FS="\t"} {a[$1]++} END {for (b in a) {print b}}' ${LIB}mapcol${INDEX14}.bowtie.txt | \
awk 'BEGIN {FS="_"; OFS="\t"} {print $1,$2,$3}' | awk '{SUM+=$3} END {if (SUM > 0) print SUM; else print "0"}')

echo $MAP14


echo "running bowtie, ${INDEX15} index, fasta input from previous unmatched, \
report perfect match, all alignments, bowtie output format"
bowtie -f --un ${LIB}unmapcol${INDEX15}.fasta -a -v 0 -p ${NSLOTS} -t ${INDEX15} ${LIB}unmapcol${INDEX14}.fasta \
> ${LIB}mapcol${INDEX15}.bowtie.txt

echo "done with mapping"

echo "head check ${LIB}mapcol${INDEX15}.bowtie.txt"
head ${LIB}mapcol${INDEX15}.bowtie.txt

echo "count how many reads were mapped for ${LIB}mapcol${INDEX15}.bowtie.txt"
MAP15=$(awk '{OFS=FS="\t"} {a[$1]++} END {for (b in a) {print b}}' ${LIB}mapcol${INDEX15}.bowtie.txt | \
awk 'BEGIN {FS="_"; OFS="\t"} {print $1,$2,$3}' | awk '{SUM+=$3} END {if (SUM > 0) print SUM; else print "0"}')

echo $MAP15



echo "running bowtie, ${INDEX16} index, fasta input from previous unmatched, \
report perfect match, all alignments, bowtie output format"
bowtie -f --un ${LIB}unmapcol${INDEX16}.fasta -a -v 0 -p ${NSLOTS} -t ${INDEX16} ${LIB}unmapcol${INDEX15}.fasta \
> ${LIB}mapcol${INDEX16}.bowtie.txt

echo "done with mapping"

echo "head check ${LIB}mapcol${INDEX16}.bowtie.txt"
head ${LIB}mapcol${INDEX16}.bowtie.txt

echo "count how many reads were mapped for ${LIB}mapcol${INDEX16}.bowtie.txt"
MAP16=$(awk '{OFS=FS="\t"} {a[$1]++} END {for (b in a) {print b}}' ${LIB}mapcol${INDEX16}.bowtie.txt | \
awk 'BEGIN {FS="_"; OFS="\t"} {print $1,$2,$3}' | awk '{SUM+=$3} END {if (SUM > 0) print SUM; else print "0"}')

echo $MAP16



echo "running bowtie, ${INDEX17} index, fasta input from previous unmatched, \
report perfect match, all alignments, bowtie output format"
bowtie -f --un ${LIB}unmapcol${INDEX17}.fasta -a -v 0 -p ${NSLOTS} -t ${INDEX17} ${LIB}unmapcol${INDEX16}.fasta \
> ${LIB}mapcol${INDEX17}.bowtie.txt

echo "done with mapping"

echo "head check ${LIB}mapcol${INDEX17}.bowtie.txt"
head ${LIB}mapcol${INDEX17}.bowtie.txt

echo "count how many reads were mapped for ${LIB}mapcol${INDEX17}.bowtie.txt"
MAP17=$(awk '{OFS=FS="\t"} {a[$1]++} END {for (b in a) {print b}}' ${LIB}mapcol${INDEX17}.bowtie.txt | \
awk 'BEGIN {FS="_"; OFS="\t"} {print $1,$2,$3}' | awk '{SUM+=$3} END {if (SUM > 0) print SUM; else print "0"}')

echo $MAP17


echo "running bowtie, ${INDEX18} index, fasta input from previous unmatched, \
report perfect match, all alignments, bowtie output format"
bowtie -f --un ${LIB}unmapcol${INDEX18}.fasta -a -v 0 -p ${NSLOTS} -t ${INDEX18} ${LIB}unmapcol${INDEX17}.fasta \
> ${LIB}mapcol${INDEX18}.bowtie.txt

echo "done with mapping"

echo "head check ${LIB}mapcol${INDEX18}.bowtie.txt"
head ${LIB}mapcol${INDEX18}.bowtie.txt

echo "count how many reads were mapped for ${LIB}mapcol${INDEX18}.bowtie.txt"
MAP18=$(awk '{OFS=FS="\t"} {a[$1]++} END {for (b in a) {print b}}' ${LIB}mapcol${INDEX18}.bowtie.txt | \
awk 'BEGIN {FS="_"; OFS="\t"} {print $1,$2,$3}' | awk '{SUM+=$3} END {if (SUM > 0) print SUM; else print "0"}')

echo $MAP18


echo "running bowtie, ${INDEX19} index, fasta input from previous unmatched, \
report perfect match, all alignments, bowtie output format"
bowtie -f --un ${LIB}unmapcol${INDEX19}.fasta -a -v 0 -p ${NSLOTS} -t ${INDEX19} ${LIB}unmapcol${INDEX18}.fasta \
> ${LIB}mapcol${INDEX19}.bowtie.txt

echo "done with mapping"

echo "head check ${LIB}mapcol${INDEX19}.bowtie.txt"
head ${LIB}mapcol${INDEX19}.bowtie.txt

echo "count how many reads were mapped for ${LIB}mapcol${INDEX19}.bowtie.txt"
MAP19=$(awk '{OFS=FS="\t"} {a[$1]++} END {for (b in a) {print b}}' ${LIB}mapcol${INDEX19}.bowtie.txt | \
awk 'BEGIN {FS="_"; OFS="\t"} {print $1,$2,$3}' | awk '{SUM+=$3} END {if (SUM > 0) print SUM; else print "0"}')

echo $MAP19


echo "running bowtie, ${INDEX20} index, fasta input from previous unmatched, \
report perfect match, all alignments, bowtie output format"
bowtie -f --un ${LIB}unmapcol${INDEX20}.fasta -a -v 0 -p ${NSLOTS} -t ${INDEX20} ${LIB}unmapcol${INDEX19}.fasta \
> ${LIB}mapcol${INDEX20}.bowtie.txt

echo "done with mapping"

echo "head check ${LIB}mapcol${INDEX20}.bowtie.txt"
head ${LIB}mapcol${INDEX20}.bowtie.txt

echo "count how many reads were mapped for ${LIB}mapcol${INDEX20}.bowtie.txt"
MAP20=$(awk '{OFS=FS="\t"} {a[$1]++} END {for (b in a) {print b}}' ${LIB}mapcol${INDEX20}.bowtie.txt | \
awk 'BEGIN {FS="_"; OFS="\t"} {print $1,$2,$3}' | awk '{SUM+=$3} END {if (SUM > 0) print SUM; else print "0"}')

echo $MAP20


echo " STEP 5 published small RNA including endo-siRNA loci"

echo "running bowtie, ${INDEX21} index, fasta input from previous unmatched, \
report perfect match, all alignments, bowtie output format"
bowtie -f --un ${LIB}unmapcol${INDEX21}.fasta -a -v 0 -p ${NSLOTS} -t ${INDEX21} ${LIB}unmapcol${INDEX20}.fasta \
> ${LIB}mapcol${INDEX21}.bowtie.txt

echo "done with mapping"

echo "head check ${LIB}mapcol${INDEX21}.bowtie.txt"
head ${LIB}mapcol${INDEX21}.bowtie.txt

echo "count how many reads were mapped for ${LIB}mapcol${INDEX21}.bowtie.txt"
MAP21=$(awk '{OFS=FS="\t"} {a[$1]++} END {for (b in a) {print b}}' ${LIB}mapcol${INDEX21}.bowtie.txt | \
awk 'BEGIN {FS="_"; OFS="\t"} {print $1,$2,$3}' | awk '{SUM+=$3} END {if (SUM > 0) print SUM; else print "0"}')

echo $MAP21


echo "running bowtie, ${INDEX22} index, fasta input from previous unmatched, \
report perfect match, all alignments, bowtie output format"
bowtie -f --un ${LIB}unmapcol${INDEX22}.fasta -a -v 0 -p ${NSLOTS} -t ${INDEX22} ${LIB}unmapcol${INDEX21}.fasta \
> ${LIB}mapcol${INDEX22}.bowtie.txt

echo "done with mapping"

echo "head check ${LIB}mapcol${INDEX22}.bowtie.txt"
head ${LIB}mapcol${INDEX22}.bowtie.txt

echo "count how many reads were mapped for ${LIB}mapcol${INDEX22}.bowtie.txt"
MAP22=$(awk '{OFS=FS="\t"} {a[$1]++} END {for (b in a) {print b}}' ${LIB}mapcol${INDEX22}.bowtie.txt | \
awk 'BEGIN {FS="_"; OFS="\t"} {print $1,$2,$3}' | awk '{SUM+=$3} END {if (SUM > 0) print SUM; else print "0"}')

echo $MAP22



echo "running bowtie, ${INDEX23} index, fasta input from previous unmatched, \
report perfect match, all alignments, bowtie output format"
bowtie -f --un ${LIB}unmapcol${INDEX23}.fasta -a -v 0 -p ${NSLOTS} -t ${INDEX23} ${LIB}unmapcol${INDEX22}.fasta \
> ${LIB}mapcol${INDEX23}.bowtie.txt

echo "done with mapping"

echo "head check ${LIB}mapcol${INDEX23}.bowtie.txt"
head ${LIB}mapcol${INDEX23}.bowtie.txt

echo "count how many reads were mapped for ${LIB}mapcol${INDEX23}.bowtie.txt"
MAP23=$(awk '{OFS=FS="\t"} {a[$1]++} END {for (b in a) {print b}}' ${LIB}mapcol${INDEX23}.bowtie.txt | \
awk 'BEGIN {FS="_"; OFS="\t"} {print $1,$2,$3}' | awk '{SUM+=$3} END {if (SUM > 0) print SUM; else print "0"}')

echo $MAP23



echo "running bowtie, ${INDEX24} index, fasta input from previous unmatched, \
report perfect match, all alignments, bowtie output format"
bowtie -f --un ${LIB}unmapcol${INDEX24}.fasta -a -v 0 -p ${NSLOTS} -t ${INDEX24} ${LIB}unmapcol${INDEX23}.fasta \
> ${LIB}mapcol${INDEX24}.bowtie.txt

echo "done with mapping"

echo "head check ${LIB}mapcol${INDEX24}.bowtie.txt"
head ${LIB}mapcol${INDEX24}.bowtie.txt

echo "count how many reads were mapped for ${LIB}mapcol${INDEX24}.bowtie.txt"
MAP24=$(awk '{OFS=FS="\t"} {a[$1]++} END {for (b in a) {print b}}' ${LIB}mapcol${INDEX24}.bowtie.txt | \
awk 'BEGIN {FS="_"; OFS="\t"} {print $1,$2,$3}' | awk '{SUM+=$3} END {if (SUM > 0) print SUM; else print "0"}')

echo $MAP24


echo "running bowtie, ${INDEX25} index, fasta input from previous unmatched, \
report perfect match, all alignments, bowtie output format"
bowtie -f --un ${LIB}unmapcol${INDEX25}.fasta -a -v 0 -p ${NSLOTS} -t ${INDEX25} ${LIB}unmapcol${INDEX24}.fasta \
> ${LIB}mapcol${INDEX25}.bowtie.txt

echo "done with mapping"

echo "head check ${LIB}mapcol${INDEX25}.bowtie.txt"
head ${LIB}mapcol${INDEX25}.bowtie.txt

echo "count how many reads were mapped for ${LIB}mapcol${INDEX25}.bowtie.txt"
MAP25=$(awk '{OFS=FS="\t"} {a[$1]++} END {for (b in a) {print b}}' ${LIB}mapcol${INDEX25}.bowtie.txt | \
awk 'BEGIN {FS="_"; OFS="\t"} {print $1,$2,$3}' | awk '{SUM+=$3} END {if (SUM > 0) print SUM; else print "0"}')

echo $MAP25


echo " STEP 6 exons overlap mRNA and non-mRNA, CISNAT"

echo "running bowtie, ${INDEX26} index, fasta input from previous unmatched, \
report perfect match, all alignments, bowtie output format"
bowtie -f --un ${LIB}unmapcol${INDEX26}.fasta -a -v 0 -p ${NSLOTS} -t ${INDEX26} ${LIB}unmapcol${INDEX25}.fasta \
> ${LIB}mapcol${INDEX26}.bowtie.txt

echo "done with mapping"

echo "head check ${LIB}mapcol${INDEX26}.bowtie.txt"
head ${LIB}mapcol${INDEX26}.bowtie.txt

echo "count how many reads were mapped for ${LIB}mapcol${INDEX26}.bowtie.txt"
MAP26=$(awk '{OFS=FS="\t"} {a[$1]++} END {for (b in a) {print b}}' ${LIB}mapcol${INDEX26}.bowtie.txt | \
awk 'BEGIN {FS="_"; OFS="\t"} {print $1,$2,$3}' | awk '{SUM+=$3} END {if (SUM > 0) print SUM; else print "0"}')

echo $MAP26



echo "running bowtie, ${INDEX27} index, fasta input from previous unmatched, \
report perfect match, all alignments, bowtie output format"
bowtie -f --un ${LIB}unmapcol${INDEX27}.fasta -a -v 0 -p ${NSLOTS} -t ${INDEX27} ${LIB}unmapcol${INDEX26}.fasta \
> ${LIB}mapcol${INDEX27}.bowtie.txt

echo "done with mapping"

echo "head check ${LIB}mapcol${INDEX27}.bowtie.txt"
head ${LIB}mapcol${INDEX27}.bowtie.txt

echo "count how many reads were mapped for ${LIB}mapcol${INDEX27}.bowtie.txt"
MAP27=$(awk '{OFS=FS="\t"} {a[$1]++} END {for (b in a) {print b}}' ${LIB}mapcol${INDEX27}.bowtie.txt | \
awk 'BEGIN {FS="_"; OFS="\t"} {print $1,$2,$3}' | awk '{SUM+=$3} END {if (SUM > 0) print SUM; else print "0"}')

echo $MAP27


echo "running bowtie, ${INDEX28} index, fasta input from previous unmatched, \
report perfect match, all alignments, bowtie output format"
bowtie -f --un ${LIB}unmapcol${INDEX28}.fasta -a -v 0 -p ${NSLOTS} -t ${INDEX28} ${LIB}unmapcol${INDEX27}.fasta \
> ${LIB}mapcol${INDEX28}.bowtie.txt

echo "done with mapping"

echo "head check ${LIB}mapcol${INDEX28}.bowtie.txt"
head ${LIB}mapcol${INDEX28}.bowtie.txt

echo "count how many reads were mapped for ${LIB}mapcol${INDEX28}.bowtie.txt"
MAP28=$(awk '{OFS=FS="\t"} {a[$1]++} END {for (b in a) {print b}}' ${LIB}mapcol${INDEX28}.bowtie.txt | \
awk 'BEGIN {FS="_"; OFS="\t"} {print $1,$2,$3}' | awk '{SUM+=$3} END {if (SUM > 0) print SUM; else print "0"}')

echo $MAP28


echo "running bowtie, ${INDEX29} index, fasta input from previous unmatched, \
report perfect match, all alignments, bowtie output format"
bowtie -f --un ${LIB}unmapcol${INDEX29}.fasta -a -v 0 -p ${NSLOTS} -t ${INDEX29} ${LIB}unmapcol${INDEX28}.fasta \
> ${LIB}mapcol${INDEX29}.bowtie.txt

echo "done with mapping"

echo "head check ${LIB}mapcol${INDEX29}.bowtie.txt"
head ${LIB}mapcol${INDEX29}.bowtie.txt

echo "count how many reads were mapped for ${LIB}mapcol${INDEX29}.bowtie.txt"
MAP29=$(awk '{OFS=FS="\t"} {a[$1]++} END {for (b in a) {print b}}' ${LIB}mapcol${INDEX29}.bowtie.txt | \
awk 'BEGIN {FS="_"; OFS="\t"} {print $1,$2,$3}' | awk '{SUM+=$3} END {if (SUM > 0) print SUM; else print "0"}')

echo $MAP29


echo "running bowtie, ${INDEX30} index, fasta input from previous unmatched, \
report perfect match, all alignments, bowtie output format"
bowtie -f --un ${LIB}unmapcol${INDEX30}.fasta -a -v 0 -p ${NSLOTS} -t ${INDEX30} ${LIB}unmapcol${INDEX29}.fasta \
> ${LIB}mapcol${INDEX30}.bowtie.txt

echo "done with mapping"

echo "head check ${LIB}mapcol${INDEX30}.bowtie.txt"
head ${LIB}mapcol${INDEX30}.bowtie.txt

echo "count how many reads were mapped for ${LIB}mapcol${INDEX30}.bowtie.txt"
MAP30=$(awk '{OFS=FS="\t"} {a[$1]++} END {for (b in a) {print b}}' ${LIB}mapcol${INDEX30}.bowtie.txt | \
awk 'BEGIN {FS="_"; OFS="\t"} {print $1,$2,$3}' | awk '{SUM+=$3} END {if (SUM > 0) print SUM; else print "0"}')

echo $MAP30


echo "count how many unmapped reads in ${LIB}unmapcol${INDEX30}.fasta"
UNMAP30=$(fasta_formatter -i ${LIB}unmapcol${INDEX30}.fasta -t | \
awk 'BEGIN {FS="_"; OFS="\t"} {print $1,$2,$3}' | awk '{SUM+=$3} END {if (SUM > 0) print SUM; else print "0"}')

echo $UNMAP30


# Remap the unmapped reads after INDEX30 to dm3excludeUext to find unique and multiple mappers

echo "running bowtie, ${MAPINDEX} index, fasta input from previous unmatched, \
report perfect match, all alignments, bowtie output format"
bowtie -f --un ${LIB}unmapcol${INDEX31}.fasta -a -v 0 -p ${NSLOTS} -t ${MAPINDEX} ${LIB}unmapcol${INDEX30}.fasta \
> ${LIB}mapcol${INDEX31}.bowtie.txt

echo "done with mapping"

echo "head check ${LIB}mapcol${INDEX31}.bowtie.txt"
head ${LIB}mapcol${INDEX31}.bowtie.txt

echo "count how many reads were mapped for ${LIB}mapcol${INDEX31}.bowtie.txt"
MAP31=$(awk '{OFS=FS="\t"} {a[$1]++} END {for (b in a) {print b}}' ${LIB}mapcol${INDEX31}.bowtie.txt | \
awk 'BEGIN {FS="_"; OFS="\t"} {print $1,$2,$3}' | awk '{SUM+=$3} END {if (SUM > 0) print SUM; else print "0"}')

echo $MAP31

echo "count how many unmapped reads in ${LIB}unmapcol${INDEX31}.fasta"
UNMAP31=$(fasta_formatter -i ${LIB}unmapcol${INDEX31}.fasta -t | \
awk 'BEGIN {FS="_"; OFS="\t"} {print $1,$2,$3}' | awk '{SUM+=$3} END {if (SUM > 0) print SUM; else print "0"}')

echo $UNMAP31



# Fix bowtie mapping file
echo "append fastaCollapseID to strand orientation"
echo "add 1 to the alignment column in bowtie standard output format"
echo "do regular sort on first column"

awk 'BEGIN {OFS=FS="\t"} {print $1,$2,$3,$4,$5,$6,$7+1}' ${LIB}mapcol${INDEX31}.bowtie.txt > ${LIB}${INDEX31}temp01.txt
sort -k1,1 ${LIB}${INDEX31}temp01.txt > ${LIB}${INDEX31}temp1.txt

echo "list uniques for fastaCollapseID and strand by removing duplicates, leaving top row"
awk '!x[$1,$2]++' ${LIB}${INDEX31}temp1.txt > ${LIB}${INDEX31}temp2.txt

echo "sum alignmentColumn+1 if have same fastaCollapseID"
awk '{OFS=FS="\t"} {a[$1]+=$7} END {for (i in a) print i,a[i]}' ${LIB}${INDEX31}temp2.txt > ${LIB}${INDEX31}temp02.txt
sort -k1,1 ${LIB}${INDEX31}temp02.txt > ${LIB}${INDEX31}temp3.txt  #This is the correct genomic hits count.

echo "convert mapping fasta file to tab format"
fasta_formatter -i ${LIB}unmapcol${INDEX30}.fasta -t > ${LIB}${INDEX31}temp03.txt
sort -k1,1 ${LIB}${INDEX31}temp03.txt > ${LIB}${INDEX31}clip${LENGTHFILTERFASTA}${FXARTIFACTFILTER}colidtab.txt

echo "do join to match mapped id with mapping fasta file"
join ${LIB}${INDEX31}temp3.txt ${LIB}${INDEX31}clip${LENGTHFILTERFASTA}${FXARTIFACTFILTER}colidtab.txt > ${LIB}${INDEX31}temp4.txt

echo "do join to match bowtie standard output alignment+1"
join ${LIB}${INDEX31}temp4.txt ${LIB}${INDEX31}temp1.txt > ${LIB}mapcol${INDEX31}.bowtie.fix.txt

echo "head check ${LIB}mapcol${INDEX31}.bowtie.fix.txt"
head ${LIB}mapcol${INDEX31}.bowtie.fix.txt


# Count unique and multiple mappers
echo "separate lists into multiple mappers and unique mappers"
awk '$2 == 1' ${LIB}mapcol${INDEX31}.bowtie.fix.txt > uniqmap${LIB}mapcol${INDEX31}.bowtie.fix.txt
awk '$2 > 1' ${LIB}mapcol${INDEX31}.bowtie.fix.txt > multimap${LIB}mapcol${INDEX31}.bowtie.fix.txt

echo "count reads that are unique mappers"
UNIMAPREMAP=$(awk '!x[$1]++' uniqmap${LIB}mapcol${INDEX31}.bowtie.fix.txt | \
awk 'BEGIN {OFS="\t"; FS="_| "} {print $3}' | \
awk '{SUM+=$1} END {if (SUM > 0) print SUM; else print "0"}')

echo $UNIMAPREMAP

echo "count reads that are multiple mappers"
MULMAPREMAP=$(awk '!x[$1]++' multimap${LIB}mapcol${INDEX31}.bowtie.fix.txt | \
awk 'BEGIN {OFS="\t"; FS="_| "} {print $3}' | \
awk '{SUM+=$1} END {if (SUM > 0) print SUM; else print "0"}')

echo $MULMAPREMAP

echo "count sum of unique and multiple mappers"
SUMREADMAPREMAPPERS=$(($UNIMAPREMAP+$MULMAPREMAP))

echo $SUMREADMAPREMAPPERS

echo "complete list"
echo $UNIMAPREMAP
echo $MULMAPREMAP
echo $SUMREADMAPREMAPPERS



echo "Complete list"
echo $LIB
echo $MAP1
echo $MAP2
echo $MAP3
echo $MAP4
echo $MAP5
echo $MAP6
echo $MAP7
echo $MAP8
echo $MAP9
echo $MAP10
echo $MAP11
echo $MAP12
echo $MAP13
echo $MAP14
echo $MAP15
echo $MAP16
echo $MAP17
echo $MAP18
echo $MAP19
echo $MAP20
echo $MAP21
echo $MAP22
echo $MAP23
echo $MAP24
echo $MAP25
echo $MAP26
echo $MAP27
echo $MAP28
echo $MAP29
echo $MAP30
echo $UNMAP30
echo $MAPGENOME
echo $UNMAPGENOME
echo $SUMREADGENOME
echo $PERGENOME
echo $UNIMAP
echo $MULMAP
echo $SUMREADMAPPERS
echo $UNIMAPREMAP
echo $MULMAPREMAP
echo $SUMREADMAPREMAPPERS


echo "rename ${LIB}unmapcol${INDEX30}.fasta to ${LIB}mapcol${INDEX31}.fasta"
mv ${LIB}unmapcol${INDEX30}.fasta ${LIB}mapcol${INDEX31}.fasta


echo "remove unneeded temp files"
rm ${LIB}temp01.txt
rm ${LIB}temp02.txt
rm ${LIB}temp03.txt
rm ${LIB}temp1.txt
rm ${LIB}temp2.txt
rm ${LIB}temp3.txt
rm ${LIB}temp4.txt
rm ${LIB}clip${LENGTHFILTERFASTA}${FXARTIFACTFILTER}colidtab.txt

rm ${LIB}mapcol${INDEXBASE}.bowtie.fix.txt
rm uniqmap${LIB}mapcol${INDEXBASE}.bowtie.fix.txt
rm multimap${LIB}mapcol${INDEXBASE}.bowtie.fix.txt

rm ${LIB}mapcol${INDEXBASE}.tab

rm ${LIB}${INDEX31}temp01.txt
rm ${LIB}${INDEX31}temp02.txt
rm ${LIB}${INDEX31}temp03.txt
rm ${LIB}${INDEX31}temp1.txt
rm ${LIB}${INDEX31}temp2.txt
rm ${LIB}${INDEX31}temp3.txt
rm ${LIB}${INDEX31}temp4.txt
rm ${LIB}${INDEX31}clip${LENGTHFILTERFASTA}${FXARTIFACTFILTER}colidtab.txt

rm ${LIB}mapcol${INDEX31}.bowtie.fix.txt
rm uniqmap${LIB}mapcol${INDEX31}.bowtie.fix.txt
rm multimap${LIB}mapcol${INDEX31}.bowtie.fix.txt

rm ${LIB}mapcol*.bowtie.txt

rm ${LIB}unmapcol*.fasta


echo "Done with ${LIB}"
date
