#!/bin/bash
#$ -cwd
#$ -V

## Join standard error and standard output into 1 file output
#$ -j y

# Time-stamp: <2016-11-23 liling>

echo "source </etc/profile.d/modules.sh>"
source /etc/profile.d/modules.sh

## Load perl module
module load perl-5.14

## Run this in the same directory with the fastq files and the list.txt.
date

echo "check version of FastQC"
perl /data/OkamuraLab/local/packagecode/FastQC/fastqc --version

preprocessfastqc.sh list.txt

date
