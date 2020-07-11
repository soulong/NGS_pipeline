
# Edit by: Hao He
# Data: 2020-01-14
# Version: 1.0

This pipeline primary was used for standard RNA sequencing analysis (including fastQC, fastq trimming, alignment(pseudot-alignment, standard alignment), and multiQC metrics).

Tools implemented in this pipeline:
fastqc (for fastq qc)
fastp (for trimming)
bowtie2 (for alignment)
STAR (for alignment)
salmon (for pseudo-alignment)
samtools (for SAM, BAM processes)
deeptools (for SAM, BAM processes)
featurecounts (for BAM quantifications)
multiqc (for multi qc metrics integration)


The prioir option to run this pipeline is though docker.
Example:

# pull docker images
docker pull antiguos/ngstools

# configuration
# edit bin/comfig.yaml, according to use definitions
# if index files were not matched to tools version, it's needed to re-index using raw fasta and gtf files

# run docker
# --volume indicates the data directory and index file directory (which can also be edited in bin/config.yaml)
docker run -it --rm --volume project_dir:/data --volume index_dir:/index antiguos/ngstools


