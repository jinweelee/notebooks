#!/bin/sh -login
#$ -cwd
#$ -V

## Join standard error and standard output into 1 file output
#$ -j y
# Time-stamp: <2016-11-23 liling>

## load perl module
module load perl-5.14

##command
date

## fetchOne - given one fastq file name, do fastQC. 
## Run and store results in current working directory. 

fetchOne() {
    FILENAME=$1
    perl /data/OkamuraLab/local/packagecode/FastQC/fastqc ${FILENAME}
    #perl /data/apps/FastQC/fastqc ${FILENAME}
}

# Do for many in fastqfilename.txt

for ARG in "$@"
do
    if [ -s "${ARG}" ]; then
        for D in `grep -v "^#" "${ARG}"`
        do
            fetchOne ${D}
        done
        exit 0
    else
        fetchOne ${ARG}
    fi
done

date
