#!/bin/bash

optspec="h:f:s:i:"

func() {
    echo "Usage:"
    echo "analysis.sh [-f config_file] [-s sample_name] [-i input_fastq_file]"
    exit -1
}

while getopts "$optspec" option; do
	case "${option}" in
		f) config_file=${OPTARG};;
		s) id=${OPTARG};;
		i) fastq=${OPTARG};;
		h) func;;
		?) func;;
	esac
done

if [ -r $config_file ]
then
	source "$config_file"
else
	echo "The config file specified: $config_file is not present."
	exit 65
fi

if [ $readlength -gt 65 ]
then
	readcut=50
	lengthcut=200
else
	readcut=30
	lengthcut=$readlength
fi

if [ ! -d $id ]
then
	mkdir $id
fi

if [ -d $id/viral_bam ]
then
	rm -r $id/viral_bam
fi

if [ -d $id/bed ]
then
	rm -r $id/bed
fi

fastp -i $fastq -y -g -o $id/$id.clean.fastq -l $readcut -w $core -h $id/$id.html

hisat2 --dta -q -x $humdb -U $id/$id.clean.fastq --sensitive -p 32 -S $id/$id.sam --un $id/$id.fastq

hisat2 --dta -q -x $virusdb -U $id/$id.fastq --bowtie2-dp 1 -k 10 --score-min L,0,-0.4 -p 32 -S $id/$id.sam

samtools view -@ $core -bS $id/$id.sam > $id/$id.bam
samtools sort -@ $core $id/$id.bam -o $id/$id.sort.bam
samtools index -@ $core $id/$id.sort.bam
samtools idxstats -@ $core $id/$id.sort.bam > $id/count_ref.txt

mkdir $id/viral_bam
mkdir $id/bed
mkdir $id/tmp
mkdir $id/single_count

python qc.py -s $id -l $lengthcut -f $annotatefile

python count_by_viral.py -s $id

python filter.py -s $id -m $method -o $outpath

if [ $method == "scRNAseq" ]
then
	python count_by_umi.py -s $id -o $outpath
	echo "done!"
fi
if [ $method == "bulkRNAseq" ]
then
	echo "done!"
fi
	
rm -r $id/tmp*
rm -r $id/single_count

rm -r $id/viral_bam

rm $id/*.fastq
rm $id/*sam
rm $id/*bam
