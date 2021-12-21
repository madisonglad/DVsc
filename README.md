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
UMI-tools  
# Create reference database
## host reference database
The first step consists in creating a hisat2 index that for host reference genomes. To do so first download the human reference genomes (GRCh38) and human mitochondrion genome.
```
mkdir /path/to/index
hisat2-build /path/to/human_reference_genomes,/path/to/human_mitochondrion_genome homo_mito
```
## virus reference database
The first step consists in creating a STAR index that include both host and virus reference genomes. To do so first download the ViruSite genome reference database. Host genome has also to be downloaded from the ensembl website. 
