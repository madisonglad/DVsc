# DVsc
DVsc is a computational software to detect and identify viruses from single-cell RNA-sequencing (scRNA-seq) or bulk RNA-seq raw data. 
This tool was tested on various scRNA-seq datasets and bulk RNA-seq datasets derived from infected human tissues and cell lines. 
DVsc was tested on a CentOS 7 cluster. 
## Dependency
The following publicly available tools should be installed to run DVsc:      
fastp https://github.com/OpenGene/fastp  
hisat2 http://daehwankimlab.github.io/hisat2/  
samtools http://www.htslib.org/   
python3 https://www.python.org/downloads/  
     --python dependecies: defaultdict, math, optparse, os, pandas, subprocess     
UMI-tools https://github.com/CGATOxford/UMI-tools   
## Create reference database
### host reference database
The first step consists in creating a hisat2 index that for host reference genomes. To do so first download the human reference genomes (GRCh38) and human mitochondrion genome from http://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz and https://www.ncbi.nlm.nih.gov/nuccore/251831106.  
To build hisat2 index, run
```
hisat2-build human_reference_genomes,human_mitochondrion_genome homodb
```
### virus reference database
The sencond step consists in creating a hisat2 index that include both host and virus reference genomes. To do so first download the virus genome reference sequences from http://www.virusite.org/index.php?nav=download and https://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=10239. You need integrate them into one file, and note that the sequence id shoud be changed as accession id, eg ```>NC_001798```.
To build hisat index, run
```
hisat2-build renamed_virus_reference_genomes,human_reference_genomes,human_mitochondrion_genome virusdb
```
## Pre-processing for scRNA-seq data
For scRNA-seq data, cellular barcode identification and UMI demultiplexing need be performed before virus detection. For droplet based techniques such as 10X and drop-seq, this can be done using UMI-tools (https://github.com/CGATOxford/UMI-tools/blob/master/doc/Single_cell_tutorial.md) with two steps
```
# step1
umi_tools whitelist --stdin Read_1.fastq --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN --log2stderr > whitelist.txt
# step2
umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN --stdin Read_1.fastq --stdout Read_extracted.fastq.gz --read2-in Read_2.fastq --read2-out=Read_extracted.fastq.gz --whitelist=whitelist.txt
```
```--bc-pattern``` should be changed depend on the sequence technique 
## Virus detection
Before run the DVsc for virus detection, you need prepare a config file, which contain all parameters. An exemple config file is provided in the Github. Then you can start the analysis by:
```
sh analysis.sh -f configfile -s sampleid -i fastqfile
```

