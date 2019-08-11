#!/bin/bash
# Bash script for trimming raw reads with Trimmomatic
# and the quality assessment of reads using FASTQC.
# Usage (call from root directory):
# scripts/trim.sh -s|-samples=[samplelist.txt] -t=[Number of Threads]
#
# files are expected to be organised as:
# root/data/
# root/data/sample1/Raw/sample1.R1.fastq.gz
# root/data/sample1/Raw/sample1.R2.fastq.gz
# root/data/sample2/Raw/sample2.fastq.gz
# etc.
# root/scripts/
# root/scripts/trim.sh
# samplelist has the name of each sample in one line.

for i in "$@"
do
case $i in
    -s=*|--samples=*)
    SAMPLELIST="${i#*=}"
    shift # past argument=value
    ;;
    -t=*|--threads=*)
    THREADS="${i#*=}"
    shift # past argument=value
    ;;
esac
done

cd data/

while read i; do
cd "$i"/Raw

### Trimmomatic ###
java -jar ~/Programs/bioinformatics/trimmomatic-0.39/dist/jar/trimmomatic-0.39.jar PE -threads ${THREADS} -phred33 *R1* *R2* ../"$i".FPs.fastq.gz ../"$i".FUs.fastq.gz  ../"$i".RPs.fastq.gz ../"$i".RUs.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:30;
echo "First sample done" "$i"
cd ..;

### FastQC ###
fastqc "$i".FPs.fastq.gz;
fastqc "$i".RPs.fastq.gz;

### unzip fastQC file ###
unzip "$i".RPs_fastqc.zip;
unzip "$i".FPs_fastqc.zip;

### Get data from FP fastqc file ###
cd "$i".FPs_fastqc/; cat fastqc_data.txt | \
 grep -E 'Filename|Total Sequences|Adapter Content|Total Deduplicated' | \
 awk '{print $NF}' | paste -sd '\t'| sed 's/.fastq//g' | \
 sed 's/FPs//' >> ../../fastqc_F.txt ; cd ..;
echo "FastQC forward" "of Sample  " "$i" "done" 
cd "$i".RPs_fastqc/; cat fastqc_data.txt | \
 grep -E 'Filename|Total Sequences|Adapter Content|Total Deduplicated' | \
 awk '{print $NF}' | paste -sd '\t'| sed 's/.fastq//g' | \
 sed 's/RPs//' >> ../../fastqc_R.txt ; cd ../;
echo "FastQC reverse" "of Sample  " "$i" "done" 
cd ..;

done < ${SAMPLELIST}
