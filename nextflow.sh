#! /bin/bash
# ubuntu_24.04_LTS
# last edited at 2025-11-13 by Hao He

#################################  setup  #################################
conda activate ngs

# download GEO
fastq-dl --group-by-sample --ignore --outdir fastq_4c --cpus 16 --force -a SRP108499 

# Hi-C
nextflow run nf-core/hic \
   -profile docker \
   --input samplesheet_hic.csv \
   --digestion dpnii \
   --fasta /mnt/f/index/hs/v49/GRCh38.primary_assembly.genome.fa \
   --chromosome_size /mnt/f/index/hs/v49/GRCh38.primary_assembly.genome.fa.gz.fai \
   --bwt2_index /mnt/f/index/hs/v49/bowtie2/ \
   --outdir results \
   --max_cpus 12 \
   --max_memory 20.GB \
   -resume

