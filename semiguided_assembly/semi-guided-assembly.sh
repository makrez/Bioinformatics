#!/bin/bash

####-----------------------####
#   Semi-guided Assembly      #
####-----------------------####

# This scripts assumes the following folder sturcture and data:
# ./samplelist.txt
# ./sample1/sample1.FPs.fastq.gz
# ./sample1/sample1.RPs.fastq.gz
# ./sample2/sample2.FPs.fastq.gz
# ./sample2/sample2.RPs.fastq.gz
# ./reference.fasta
#
# list of samples in samplelist.txt (e.g. sample1 \n sample2 \n ...)

# Program Paths (need to be adjusted). samtools, bwa, spades etc. are in PATH.

BWA=bwa
PICARD=~/Programs/bioinformatics/picard.jar

for i in "$@"
do
case $i in
    -s=*|--samples=*)
    SAMPLELIST="${i#*=}"
    shift # past argument=value
    ;;
    -r=*|--ref=*)
    REFERENCE="${i#*=}"
    shift # past argument=value
    ;;
    -o=*|--output_name=*)
    OUTNAME="${i#*=}"
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

## Status of mapping
cd $i; echo $i >> ../"status_mapping${OUTNAME}.txt";

#### Mapping reads to reference ####
echo "1 Mapping starting"
${BWA} mem -t ${THREADS} -M ../${REFERENCE} "$i".FPs.fastq.gz "$i".RPs.fastq.gz > "$i".${OUTNAME}.sam;
echo "2 BWA mem done"
samtools view -b "$i".${OUTNAME}.sam > "$i".${OUTNAME}.bam;
echo "3 Samtools view done"
samtools sort --threads ${THREADS} "$i".${OUTNAME}.bam > $i.${OUTNAME}.sorted; echo $i >> ../stats_${OUTNAME}.txt;
echo "4 Samtools sort done"
samtools stats -d $i.${OUTNAME}.sorted | grep ^SN | cut -f 2- >> ../stats_${OUTNAME}_mapping.txt;
echo "5 Samtools stats done"
sudo rm $i.${OUTNAME}.sam;
sudo rm $i.${OUTNAME}.bam;
echo "Mapping sample finished"
cd ..; done < ../${SAMPLELIST}

# De novo assembly using mapped reads to a reference #
# Sorted bam file is already present #

while read i; do
cd "$i";
echo "denovo starting"
samtools view -b -F 4 "$i".${OUTNAME}.sorted > "$i".${OUTNAME}.mapped.bam;

java -jar $PICARD SamToFastq I="$i".${OUTNAME}.mapped.bam F="$i".${OUTNAME}.F.fastq F2="$i".${OUTNAME}.R.fastq;

### Assembly of Reads ###

spades.py --careful -t ${THREADS} -o ${OUTNAME}.assembly -1 "$i".${OUTNAME}.F.fastq -2 "$i".${OUTNAME}.R.fastq;
echo "Spades finished"
cd ${OUTNAME}.assembly; cat contigs.fasta | awk  '/^>/{if(N)exit;++N;} {print;}' | sed "s/NODE_1/$i/g" > "$i".${OUTNAME}.fasta;
cat "$i".${OUTNAME}.fasta | grep '>' | awk -v var="$i" -F '_' '{print var"\t"$6}' >> ../../${OUTNAME}_mapping_coverage.txt; cd ..;
cd ..; done < ../${SAMPLELIST}
