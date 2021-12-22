# DVsc
DVsc is a computational software to detect and identify viruses from single-cell RNA-sequencing (scRNA-seq) or bulk RNA-seq raw data. 
This tool was tested on various scRNA-seq datasets and bulk RNA-seq datasets derived from infected human tissues and cell lines. 
DVsc was tested on a CentOS 7 cluster. 
# Dependency
The following publicly available tools should be installed to run DVsc:    
fastp   
hisat2  
samtools  
python3
--依赖的相关模块()
UMI-tools  
# Create reference database
## host reference database
The first step consists in creating a hisat2 index that for host reference genomes. To do so first download the human reference genomes (GRCh38) and human mitochondrion genome from http://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz and https://www.ncbi.nlm.nih.gov/nuccore/251831106.  
To build hisat index, run
```
mkdir /path/to/index
hisat2-build /path/to/human_reference_genomes,/path/to/human_mitochondrion_genome homo_mito
```
## virus reference database
The sencond step consists in creating a hisat2 index that include both host and virus reference genomes. To do so first download the virus genome reference sequences from http://www.virusite.org/index.php?nav=download and https://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=10239. You need integrate them into one file, and rename sequence id using (.py).  
To build hisat index, run
```
mkdir /path/to/index
hisat2-build /path/to/renamed_virus_reference_genomes,/path/to/human_reference_genomes,/path/to/human_mitochondrion_genome homo_mito
```
