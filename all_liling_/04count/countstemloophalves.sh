#!/bin/bash
# Time-stamp: <2015-12-10 08:54:28 liling>

##################################################################################################
# This script does the job of counting the number of rawcount and normreadcount
# sequences that mapped to the adjusted stemloop halves bed file for dme mirbasev21plusITS1.
##################################################################################################

## exit when there is an error
set -e

## set locale as C
export LC_ALL="C"

echo "check locale setting"
locale

echo "Start"
date

echo "This is running countstemloophalves script from $PWD"

printf "\n"

echo "input file is a continuation from the countfeaturesbysensereads script"

for LIB in B133 B134 B135 B136 B137 B138 B139 B140 B141 B142 B143 B144 B145 ; do
    for INDEXBASE in dmestemloopmirbasev21plusITS1 ; do

        echo "head check ${LIB}mapcol${INDEXBASE}.bowtie.txtsense.fix.maphitdm6.rawcount.normreadcount.txt"
        echo "<identifier> <maphitdm6> <lib> <rawcount> <normreadcount> <orientation> <geneID> <start> <bowtiesequence> <bowtiequal> <bowtiehit>"
        head ${LIB}mapcol${INDEXBASE}.bowtie.txtsense.fix.maphitdm6.rawcount.normreadcount.txt

        printf "\n"

        echo "adjust to BED format with the different value counts"
        echo "do for rawcount"
        echo "CHR START STOP <:${LIB}:sequence:countvalue:mapHit:lengthnt> COUNTVALUE ORIENTATION"
        awk 'BEGIN {OFS=FS="\t"} {print $7,$8,($8+length($9)),$3":"$9":"$4":"$2":"length($9)"nt",$4,$6}' ${LIB}mapcol${INDEXBASE}.bowtie.txtsense.fix.maphitdm6.rawcount.normreadcount.txt \
            > ${LIB}mapcol${INDEXBASE}.sense.maphitdm6.rawcount.bed

        echo "use sort-bed to sort bed files"
        sort-bed --max-mem 4G --tmpdir $PWD ${LIB}mapcol${INDEXBASE}.sense.maphitdm6.rawcount.bed > ${LIB}mapcol${INDEXBASE}.sense.maphitdm6.rawcount.bedsort
        mv ${LIB}mapcol${INDEXBASE}.sense.maphitdm6.rawcount.bedsort ${LIB}mapcol${INDEXBASE}.sense.maphitdm6.rawcount.bed

        echo "head check ${LIB}mapcol${INDEXBASE}.sense.maphitdm6.rawcount.bed"
        echo "CHR START STOP <:${LIB}:sequence:countvalue:mapHit:lengthnt> RAWCOUNT ORIENTATION"
        head ${LIB}mapcol${INDEXBASE}.sense.maphitdm6.rawcount.bed
        wc -l ${LIB}mapcol${INDEXBASE}.sense.maphitdm6.rawcount.bed

        printf "\n"

        echo "do for normreadcount"
        echo "CHR START STOP <:${LIB}:sequence:countvalue:mapHit:lengthnt> COUNTVALUE ORIENTATION"
        awk 'BEGIN {OFS=FS="\t"} {print $7,$8,($8+length($9)),$3":"$9":"$5":"$2":"length($9)"nt",$5,$6}' ${LIB}mapcol${INDEXBASE}.bowtie.txtsense.fix.maphitdm6.rawcount.normreadcount.txt \
            > ${LIB}mapcol${INDEXBASE}.sense.maphitdm6.normreadcount.bed

        echo "use sort-bed to sort bed files"
        sort-bed --max-mem 4G --tmpdir $PWD ${LIB}mapcol${INDEXBASE}.sense.maphitdm6.normreadcount.bed > ${LIB}mapcol${INDEXBASE}.sense.maphitdm6.normreadcount.bedsort
        mv ${LIB}mapcol${INDEXBASE}.sense.maphitdm6.normreadcount.bedsort ${LIB}mapcol${INDEXBASE}.sense.maphitdm6.normreadcount.bed

        echo "head check ${LIB}mapcol${INDEXBASE}.sense.maphitdm6.normreadcount.bed"
        echo "CHR START STOP <:${LIB}:sequence:countvalue:mapHit:lengthnt> NORMREADCOUNT ORIENTATION"
        head ${LIB}mapcol${INDEXBASE}.sense.maphitdm6.normreadcount.bed
        wc -l ${LIB}mapcol${INDEXBASE}.sense.maphitdm6.normreadcount.bed

        printf "\n"
    done
    printf "\n"
done

echo "use sort-bed to sort bed files"
sort-bed --max-mem 4G --tmpdir $PWD dmestemloopmirbasev21plusITS1halves.bed > dmestemloopmirbasev21plusITS1halves.bedsort
mv dmestemloopmirbasev21plusITS1halves.bedsort dmestemloopmirbasev21plusITS1halves.bed

echo "head check dmestemloopmirbasev21plusITS1halves.bed"
head dmestemloopmirbasev21plusITS1halves.bed
wc -l dmestemloopmirbasev21plusITS1halves.bed

printf "\n"

for LIB in B133 B134 B135 B136 B137 B138 B139 B140 B141 B142 B143 B144 B145 ; do
    for INDEXBASE in dmestemloopmirbasev21plusITS1 ; do

        echo "count reads that mapped exactly within interval region in dmestemloopmirbasev21plusITS1halves.bed for rawcount"
        bedtools map -a dmestemloopmirbasev21plusITS1halves.bed -b ${LIB}mapcol${INDEXBASE}.sense.maphitdm6.rawcount.bed -c 5 -o sum -F 1.0 -null "0" \
            > ${LIB}.dmestemloopmirbasev21plusITS1halves.rawcount.txt

        echo "use sort-bed to sort files"
        sort-bed --max-mem 4G --tmpdir $PWD ${LIB}.dmestemloopmirbasev21plusITS1halves.rawcount.txt > ${LIB}.dmestemloopmirbasev21plusITS1halves.rawcount.txtsort
        mv ${LIB}.dmestemloopmirbasev21plusITS1halves.rawcount.txtsort ${LIB}.dmestemloopmirbasev21plusITS1halves.rawcount.txt

        echo "head check ${LIB}.dmestemloopmirbasev21plusITS1halves.rawcount.txt"
        head ${LIB}.dmestemloopmirbasev21plusITS1halves.rawcount.txt
        wc -l ${LIB}.dmestemloopmirbasev21plusITS1halves.rawcount.txt

        printf "\n"
        
        echo "print only value"
        echo "add header"
        awk 'BEGIN {OFS=FS="\t"} {print $7}' ${LIB}.dmestemloopmirbasev21plusITS1halves.rawcount.txt > ${LIB}.dmestemloopmirbasev21plusITS1halves.rawcount.value.txt
        sed -i $'1 i\\\n'$LIB'.rawcount' ${LIB}.dmestemloopmirbasev21plusITS1halves.rawcount.value.txt


        printf "\n"

        echo "count reads that mapped exactly within interval region in dmestemloopmirbasev21plusITS1halves.bed for normreadcount"
        bedtools map -a dmestemloopmirbasev21plusITS1halves.bed -b ${LIB}mapcol${INDEXBASE}.sense.maphitdm6.normreadcount.bed -c 5 -o sum -F 1.0 -null "0" \
            > ${LIB}.dmestemloopmirbasev21plusITS1halves.normreadcount.txt

        echo "use sort-bed to sort files"
        sort-bed --max-mem 4G --tmpdir $PWD ${LIB}.dmestemloopmirbasev21plusITS1halves.normreadcount.txt > ${LIB}.dmestemloopmirbasev21plusITS1halves.normreadcount.txtsort
        mv ${LIB}.dmestemloopmirbasev21plusITS1halves.normreadcount.txtsort ${LIB}.dmestemloopmirbasev21plusITS1halves.normreadcount.txt

        echo "head check ${LIB}.dmestemloopmirbasev21plusITS1halves.normreadcount.txt"
        head ${LIB}.dmestemloopmirbasev21plusITS1halves.normreadcount.txt
        wc -l ${LIB}.dmestemloopmirbasev21plusITS1halves.normreadcount.txt

        printf "\n"

        echo "print only value"
        echo "add header"
        awk 'BEGIN {OFS=FS="\t"} {print $7}' ${LIB}.dmestemloopmirbasev21plusITS1halves.normreadcount.txt > ${LIB}.dmestemloopmirbasev21plusITS1halves.normreadcount.value.txt
        sed -i $'1 i\\\n'$LIB'.normreadcount' ${LIB}.dmestemloopmirbasev21plusITS1halves.normreadcount.value.txt

    done
    printf "\n"
done
printf "\n"

echo "add header to dmestemloopmirbasev21plusITS1halves list"
sed $'1 i\\\nMIRID\tSTART\tSTOP\tSTEMLOOPHALF\tFAKESCORE\tORIENTATION' dmestemloopmirbasev21plusITS1halves.bed > dmestemloopmirbasev21plusITS1halves.header.txt

echo "head check dmestemloopmirbasev21plusITS1halves.header.txt"
head dmestemloopmirbasev21plusITS1halves.header.txt
wc -l dmestemloopmirbasev21plusITS1halves.header.txt

echo "combine files together for tables"
echo "do for rawcount"
paste *.dmestemloopmirbasev21plusITS1halves.rawcount.value.txt > dmestemloopmirbasev21plusITS1halves.rawcount.value.txttemp
paste dmestemloopmirbasev21plusITS1halves.header.txt dmestemloopmirbasev21plusITS1halves.rawcount.value.txttemp > dmestemloopmirbasev21plusITS1halves.rawcount.table.txt
rm dmestemloopmirbasev21plusITS1halves.rawcount.value.txttemp

echo "head check dmestemloopmirbasev21plusITS1halves.rawcount.table.txt"
head dmestemloopmirbasev21plusITS1halves.rawcount.table.txt

printf "\n"

echo "do for normreadcount"
paste *.dmestemloopmirbasev21plusITS1halves.normreadcount.value.txt > dmestemloopmirbasev21plusITS1halves.normreadcount.value.txttemp
paste dmestemloopmirbasev21plusITS1halves.header.txt dmestemloopmirbasev21plusITS1halves.normreadcount.value.txttemp > dmestemloopmirbasev21plusITS1halves.normreadcount.table.txt
rm dmestemloopmirbasev21plusITS1halves.normreadcount.value.txttemp

echo "head check dmestemloopmirbasev21plusITS1halves.normreadcount.table.txt"
head dmestemloopmirbasev21plusITS1halves.normreadcount.table.txt

printf "\n"

echo "Done"
date
