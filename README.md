# Streptococcus suis whole genome assembly and typing pipeline


Pipeline for whole genome assembly and analysis of Streptococccus suis. Works only for Oxford nanopore reads


Pipeline is based on this paper 
Athey, T.B.T., Teatero, S., Lacouture, S. et al. Determining Streptococcus suis serotype from short-read whole-genome sequencing data. BMC Microbiol 16, 162 (2016). https://doi.org/10.1186/s12866-016-0782-8

Requires input directory containg sub-directories with the fastq files and output directory. Outputs several intermediate files with a html report with AMR,MLST,Serotyping and virulence factors found in the sample. 

Usage
```
nextflow run main.nf --input sampleslist.csv --outdir Results_mannheimia_2 -profile docker --trim_barcodes
```
```
Parameters:

--input		Input directory containg sub-sirectories with fastq files
--out_dir	Output directory
optional
--trim_barcodes barcode and adapter trimming using porechop
```
## Dependencies
* nextflow
* docker
* wsl2
## Software and references used
* dragonflye (https://github.com/rpetit3/dragonflye)
* Database for serotyping and virulence typing was constrcuted using sequences from https://github.com/streplab/SsuisSerotyping_pipeline 
* porechop (https://github.com/rrwick/Porechop)
* minimap2 (Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34:3094-3100. doi:10.1093/bioinformatics/bty191)
* abricate (https://github.com/tseemann/abricate)
* mlst (https://github.com/tseemann/mlst,This publication made use of the PubMLST website (https://pubmlst.org/) developed by Keith Jolley (Jolley & Maiden 2010, BMC        Bioinformatics, 11:595) and sited at the University of Oxford. The development of that website was funded by the Wellcome Trust)
* rmarkdown 
* samtools and bcftools (http://www.htslib.org/)
