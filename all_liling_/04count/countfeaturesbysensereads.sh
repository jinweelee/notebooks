#!/bin/bash
# Time-stamp: <2015-12-04 18:24:46 liling>


##################################################################################################################################
# LIB = library name
# INDEXBASE = bowtie mapping index, same as in the mapcolfasta or mapcolfastafeatures script
# ASSEMBLY = UCSC assembly identifier, eg. mm10, dm6, hg19
##################################################################################################################################


## exit when there is an error
set -e

## set locale as C
export LC_ALL="C"

echo "check locale setting"
locale

echo "Start"
date

echo "This is running countfeaturesbysensereads script from $PWD"

printf "\n"

for LIB in B133 B134 B135 B136 B137 B138 B139 B140 B141 B142 B143 B144 B145 ; do

    ##################################################################################################################################
    # Find the mapHit for each identifer in dm6 from the bowtie output file.
    # Do counting fix to get the mapHits within the index itself
    # note that there might be some mapping across different categories
    # count number of times a read identifier appears in the $INPUTBOWTIEFILE
    # sort by first key

    echo "head check ${LIB}mapcoldm6.bowtie.txt"
    head ${LIB}mapcoldm6.bowtie.txt

    printf "\n"

    echo "Count number of times a read identifier appears in the ${LIB}mapcoldm6.bowtie.txt"
    echo "note that the identifier can map on either strand"
    awk 'BEGIN {OFS=FS="\t"} {a[$1]++} END {for (i in a) {print i,a[i]}}' ${LIB}mapcoldm6.bowtie.txt \
        > ${LIB}mapcoldm6.bowtie.txt.idlist.maphitdm6.txt

    echo "sort by first key"
    sort -k1,1 ${LIB}mapcoldm6.bowtie.txt.idlist.maphitdm6.txt > ${LIB}mapcoldm6.bowtie.txt.idlist.maphitdm6.txtsort

    echo "format to tab file"
    awk 'BEGIN {OFS="\t"; FS="_|\t"} {print $1"_"$2"_"$3,$4,$1,$3,sprintf("%.4f",($3/$4))}' ${LIB}mapcoldm6.bowtie.txt.idlist.maphitdm6.txtsort \
        > ${LIB}mapcoldm6.bowtie.txt.idlist.maphitdm6.txt

    rm ${LIB}mapcoldm6.bowtie.txt.idlist.maphitdm6.txtsort

    printf "\n"

    echo "head check ${LIB}mapcoldm6.bowtie.txt.idlist.maphitdm6.txt"
    echo "<IDENTIFIER> <MAPHITdm6> <LIB> <RAWCOUNT> <NORMREADCOUNT>"
    head ${LIB}mapcoldm6.bowtie.txt.idlist.maphitdm6.txt
    wc -l ${LIB}mapcoldm6.bowtie.txt.idlist.maphitdm6.txt

    printf "\n"
    
    echo "Sum the RAWCOUNT in ${LIB}mapcoldm6.bowtie.txt.idlist.maphitdm6.txt"
    awk '{SUM+=$4} END {if (SUM > 0) printf("%.4f\n", SUM); else print "0"}' ${LIB}mapcoldm6.bowtie.txt.idlist.maphitdm6.txt

    printf "\n"

done

printf "\n"


##################################################################################################################################
# Select only matches to sense strand for the features mapped files

for LIB in B133 B134 B135 B136 B137 B138 B139 B140 B141 B142 B143 B144 B145 ; do
    for INDEXBASE in spikeinset01 dmestemloopmirbasev21plusITS1 ; do

        echo "head check ${LIB}mapcol${INDEXBASE}.bowtie.txt"
        head ${LIB}mapcol${INDEXBASE}.bowtie.txt

        printf "\n"

        echo "Select only matches to sense strand"
        awk '$2 == "+"' ${LIB}mapcol${INDEXBASE}.bowtie.txt > ${LIB}mapcol${INDEXBASE}.bowtie.txtsense

        echo "head check ${LIB}mapcol${INDEXBASE}.bowtie.txtsense"
        head ${LIB}mapcol${INDEXBASE}.bowtie.txtsense

        printf "\n"

        echo "count how many reads mapped to the index in sense"
        echo "consider each identifier only once"
        echo "note that identifiers may be mapped to more than one geneID"
        LIBMAPSUMSENSE=$(awk 'BEGIN {OFS=FS="\t"} {a[$1]++} END {for (b in a) {print b}}' ${LIB}mapcol${INDEXBASE}.bowtie.txtsense | \
            awk 'BEGIN {OFS="\t"; FS="_"} {print $1,$2,$3}' | awk '{SUM+=$3} END {if (SUM > 0) print SUM; else print "0"}')

        echo $LIBMAPSUMSENSE

        printf "\n"

        echo "sort ${LIB}mapcol${INDEXBASE}.bowtie.txtsense by first key"
        sort -k1,1 ${LIB}mapcol${INDEXBASE}.bowtie.txtsense > ${LIB}mapcol${INDEXBASE}.bowtie.txtsensesort
        mv ${LIB}mapcol${INDEXBASE}.bowtie.txtsensesort ${LIB}mapcol${INDEXBASE}.bowtie.txtsense

        echo "head check ${LIB}mapcol${INDEXBASE}.bowtie.txtsense"
        head ${LIB}mapcol${INDEXBASE}.bowtie.txtsense
        
        printf "\n"

    done
    printf "\n"
done

printf "\n"

for LIB in B133 B134 B135 B136 B137 B138 B139 B140 B141 B142 B143 B144 B145 ; do
    for INDEXBASE in dmestemloopmirbasev21plusITS1 ; do
        
        echo "join ${LIB}mapcoldm6.bowtie.txt.idlist.maphitdm6.txt with ${LIB}mapcol${INDEXBASE}.bowtie.txtsense to get the maphit matched"
        join -t $'\t' ${LIB}mapcoldm6.bowtie.txt.idlist.maphitdm6.txt ${LIB}mapcol${INDEXBASE}.bowtie.txtsense \
            > ${LIB}mapcol${INDEXBASE}.bowtie.txtsense.fix.maphitdm6.rawcount.normreadcount.txt
        
        echo "head check ${LIB}mapcol${INDEXBASE}.bowtie.txtsense.fix.maphitdm6.rawcount.normreadcount.txt"
        echo "<identifier> <maphitdm6> <lib> <rawcount> <normreadcount> <orientation> <geneID> <start> <bowtiesequence> <bowtiequal> <bowtiehit>"
        head ${LIB}mapcol${INDEXBASE}.bowtie.txtsense.fix.maphitdm6.rawcount.normreadcount.txt

        echo "Calculate sum of values"
        echo "count rawcount for each geneID in ${INDEXBASE}"
        awk 'BEGIN {OFS=FS="\t"} {a[$7]+=$4} END {for (i in a) printf("%s\t%.4f\n", i,a[i])}' ${LIB}mapcol${INDEXBASE}.bowtie.txtsense.fix.maphitdm6.rawcount.normreadcount.txt \
            > ${LIB}.${INDEXBASE}.rawcount.txt
        sort -k1,1 ${LIB}.${INDEXBASE}.rawcount.txt > ${LIB}.${INDEXBASE}.rawcount.txtsort
        mv ${LIB}.${INDEXBASE}.rawcount.txtsort ${LIB}.${INDEXBASE}.rawcount.txt

        echo "head check ${LIB}.${INDEXBASE}.rawcount.txt"
        echo "<GENEID> <RAWCOUNT>"
        head ${LIB}.${INDEXBASE}.rawcount.txt
        echo "count SUM"
        awk '{SUM+=$2} END {if (SUM > 0) printf("%.4f\n", SUM); else print "0"}' ${LIB}.${INDEXBASE}.rawcount.txt
        echo "count number of lines"
        wc -l ${LIB}.${INDEXBASE}.rawcount.txt

        printf "\n"

        echo "count normreadcount for each geneID in ${INDEXBASE}"
        awk 'BEGIN {OFS=FS="\t"} {a[$7]+=$5} END {for (i in a) printf("%s\t%.4f\n", i,a[i])}' ${LIB}mapcol${INDEXBASE}.bowtie.txtsense.fix.maphitdm6.rawcount.normreadcount.txt \
            > ${LIB}.${INDEXBASE}.normreadcount.txt
        sort -k1,1 ${LIB}.${INDEXBASE}.normreadcount.txt > ${LIB}.${INDEXBASE}.normreadcount.txtsort
        mv ${LIB}.${INDEXBASE}.normreadcount.txtsort ${LIB}.${INDEXBASE}.normreadcount.txt

        echo "head check ${LIB}.${INDEXBASE}.normreadcount.txt"
        echo "<GENEID> <NORMREADCOUNT>"
        head ${LIB}.${INDEXBASE}.normreadcount.txt
        echo "count SUM"
        awk '{SUM+=$2} END {if (SUM > 0) printf("%.4f\n", SUM); else print "0"}' ${LIB}.${INDEXBASE}.normreadcount.txt
        echo "count number of lines"
        wc -l ${LIB}.${INDEXBASE}.normreadcount.txt

        printf "\n"

    done
    printf "\n"
done

printf "\n"


for LIB in B133 B134 B135 B136 B137 B138 B139 B140 B141 B142 B143 B144 B145 ; do
    for INDEXBASE in spikeinset01 ; do
        
        echo "Calculate sum of values"
        echo "note that spikeinset01 does not match to genome, so, no mapHitdm6"

        printf "\n"

        echo "adjust ${LIB}mapcol${INDEXBASE}.bowtie.txtsense file"
        awk 'BEGIN {OFS=FS="\t"} {split($1,id,"_"); print $1,id[1],id[3],$2,$3,$4,$5}' ${LIB}mapcol${INDEXBASE}.bowtie.txtsense > ${LIB}mapcol${INDEXBASE}.bowtie.txtsense.rawcount.txt

        printf "\n"

        echo "head check ${LIB}mapcol${INDEXBASE}.bowtie.txtsense.rawcount.txt"
        echo "<identifier> <lib> <rawcount> <orientation> <geneid> <start> <bowtiesequence>"
        head ${LIB}mapcol${INDEXBASE}.bowtie.txtsense.rawcount.txt

        printf "\n"

        echo "count rawcount for each geneID in ${LIB}mapcol${INDEXBASE}.bowtie.txtsense"
        awk 'BEGIN {OFS=FS="\t"} {a[$5]+=$3} END {for (i in a) printf("%s\t%.4f\n", i,a[i])}' ${LIB}mapcol${INDEXBASE}.bowtie.txtsense.rawcount.txt \
            > ${LIB}.${INDEXBASE}.rawcount.txt
        sort -k1,1 ${LIB}.${INDEXBASE}.rawcount.txt > ${LIB}.${INDEXBASE}.rawcount.txtsort
        mv ${LIB}.${INDEXBASE}.rawcount.txtsort ${LIB}.${INDEXBASE}.rawcount.txt

        echo "head check ${LIB}.${INDEXBASE}.rawcount.txt"
        echo "<GENEID> <RAWCOUNT>"
        head ${LIB}.${INDEXBASE}.rawcount.txt
        echo "count SUM"
        awk '{SUM+=$2} END {if (SUM > 0) printf("%.4f\n", SUM); else print "0"}' ${LIB}.${INDEXBASE}.rawcount.txt
        echo "count number of lines"
        wc -l ${LIB}.${INDEXBASE}.rawcount.txt

        printf "\n"

   done
    printf "\n"
done

echo "perform sequential join with the full index list"
echo "add 0 for empty entries"

printf "\n"

echo "format index list"

for INDEXBASE in spikeinset01 dmestemloopmirbasev21plusITS1 ; do

    printf "\n"

    echo "sort full index list with same sort key"
    sort -k1,1 ${INDEXBASE}.list.txt > ${INDEXBASE}.list.txtsort
    mv ${INDEXBASE}.list.txtsort ${INDEXBASE}.list.txt

    echo "head check ${INDEXBASE}.list.txt"
    echo "<GENEID> <ACCESSION> <LENGTH>"
    head ${INDEXBASE}.list.txt
    wc -l ${INDEXBASE}.list.txt

    echo "add header for the sorted indexlist"
    sed $'1 i\\\nGENEID\tACCESSION\tLENGTH' ${INDEXBASE}.list.txt > ${INDEXBASE}.list.txt.header 

    echo "head check ${INDEXBASE}.list.txt.header"
    head ${INDEXBASE}.list.txt.header
    wc -l ${INDEXBASE}.list.txt.header

    printf "\n"

done

printf "\n"

echo "do for count values"
for LIB in B133 B134 B135 B136 B137 B138 B139 B140 B141 B142 B143 B144 B145 ; do
    for INDEXBASE in dmestemloopmirbasev21plusITS1 ; do
        for COUNTTYPE in rawcount normreadcount ; do
            
            printf "\n"

            echo "do join for ${LIB} ${INDEXBASE} ${COUNTTYPE}"
            echo "<GENEID> <COUNTVALUE>"
            join -a 1 -t $'\t' -o 0,2.2 -e "0" ${INDEXBASE}.list.txt ${LIB}.${INDEXBASE}.${COUNTTYPE}.txt > ${LIB}.${INDEXBASE}.${COUNTTYPE}.join.txt

            head ${LIB}.${INDEXBASE}.${COUNTTYPE}.join.txt

            echo "print only value"
            echo "add header"
            awk 'BEGIN {OFS=FS="\t"} {print $2}' ${LIB}.${INDEXBASE}.${COUNTTYPE}.join.txt > ${LIB}.${INDEXBASE}.${COUNTTYPE}.value.txt
            sed -i $'1 i\\\n'$LIB'.'$INDEXBASE'.'$COUNTTYPE'' ${LIB}.${INDEXBASE}.${COUNTTYPE}.value.txt

            head ${LIB}.${INDEXBASE}.${COUNTTYPE}.value.txt
            
            printf "\n"
        done
    done
done

printf "\n"

echo "do for count values"
for LIB in B133 B134 B135 B136 B137 B138 B139 B140 B141 B142 B143 B144 B145 ; do
    for INDEXBASE in spikeinset01 ; do
        for COUNTTYPE in rawcount ; do
            
            printf "\n"

            echo "do join for ${LIB} ${INDEXBASE} ${COUNTTYPE}"
            echo "<GENEID> <COUNTVALUE>"
            join -a 1 -t $'\t' -o 0,2.2 -e "0" ${INDEXBASE}.list.txt ${LIB}.${INDEXBASE}.${COUNTTYPE}.txt > ${LIB}.${INDEXBASE}.${COUNTTYPE}.join.txt

            head ${LIB}.${INDEXBASE}.${COUNTTYPE}.join.txt

            echo "print only value"
            echo "add header"
            awk 'BEGIN {OFS=FS="\t"} {print $2}' ${LIB}.${INDEXBASE}.${COUNTTYPE}.join.txt > ${LIB}.${INDEXBASE}.${COUNTTYPE}.value.txt
            sed -i $'1 i\\\n'$LIB'.'$INDEXBASE'.'$COUNTTYPE'' ${LIB}.${INDEXBASE}.${COUNTTYPE}.value.txt

            head ${LIB}.${INDEXBASE}.${COUNTTYPE}.value.txt
            
            printf "\n"
        done
    done
done

printf "\n"

for INDEXBASE in dmestemloopmirbasev21plusITS1 ; do

    echo "combine together files for tables ${INDEXBASE}"
    echo "do for rawcount"
    paste *.${INDEXBASE}.rawcount.value.txt > ${INDEXBASE}.rawcount.value.txttemp1
    paste ${INDEXBASE}.list.txt.header ${INDEXBASE}.rawcount.value.txttemp1 > ${INDEXBASE}.rawcount.table.txt
    
    printf "\n"

    echo "do for normreadcount"
    paste *.${INDEXBASE}.normreadcount.value.txt > ${INDEXBASE}.normreadcount.value.txttemp1
    paste ${INDEXBASE}.list.txt.header ${INDEXBASE}.normreadcount.value.txttemp1 > ${INDEXBASE}.normreadcount.table.txt

    rm ${INDEXBASE}.rawcount.value.txttemp1
    rm ${INDEXBASE}.normreadcount.value.txttemp1

    printf "\n"
done

printf "\n"

for INDEXBASE in spikeinset01 ; do

    echo "combine together files for tables ${INDEXBASE}"
    echo "do for rawcount"
    paste *.${INDEXBASE}.rawcount.value.txt > ${INDEXBASE}.rawcount.value.txttemp1
    paste ${INDEXBASE}.list.txt.header ${INDEXBASE}.rawcount.value.txttemp1 > ${INDEXBASE}.rawcount.table.txt
    
    printf "\n"

    rm ${INDEXBASE}.rawcount.value.txttemp1

    printf "\n"
done

printf "\n"

echo "do cleanup"
for LIB in B133 B134 B135 B136 B137 B138 B139 B140 B141 B142 B143 B144 B145 ; do
    for INDEXBASE in dmestemloopmirbasev21plusITS1 ; do
        for COUNTTYPE in rawcount normreadcount ; do

            echo "remove unneeded files"
            rm ${LIB}.${INDEXBASE}.${COUNTTYPE}.txt
            rm ${LIB}.${INDEXBASE}.${COUNTTYPE}.join.txt
            rm ${LIB}.${INDEXBASE}.${COUNTTYPE}.value.txt

            printf "\n"
        done
    done
done


echo "do cleanup"
for LIB in B133 B134 B135 B136 B137 B138 B139 B140 B141 B142 B143 B144 B145 ; do
    for INDEXBASE in spikeinset01 ; do
        for COUNTTYPE in rawcount ; do

            echo "remove unneeded files"
            rm ${LIB}.${INDEXBASE}.${COUNTTYPE}.txt
            rm ${LIB}.${INDEXBASE}.${COUNTTYPE}.join.txt
            rm ${LIB}.${INDEXBASE}.${COUNTTYPE}.value.txt

            printf "\n"
        done
    done
done

printf "\n"

echo "This done running countfeaturesbysensereads script from $PWD"

echo "Done"
date
