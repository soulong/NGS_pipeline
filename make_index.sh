#! /bin/bash
# ubuntu_26.04_LTS
# last edited at 2026-07-08 by Hao He

# conda activate ngs
# conda env export --no-builds | grep -v "^prefix: " > ngs.yml
# conda env create -f ngs.yml

#################################  setup  #################################

## T2T-CHM13v2.0 (https://github.com/marbl/CHM13), neeed to further processed by scripts before indexing
# fasta_link=https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
# gtf_link=https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_RefSeq_Liftoff_v5.3.gff.gz

## hs, chromosome: chr1, ..., chrX, chrY, chrM, GLxxx, KIxxx
# fasta_link=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.primary_assembly.genome.fa.gz
# gtf_link=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.basic.annotation.gtf.gz
# cdna_link=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.transcripts.fa.gz

## mm, chromosome: chr1, ..., chrX, chrY, chrM, GLxxx, JHxxx, MUxxx
# fasta_link=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/GRCm39.primary_assembly.genome.fa.gz
# cdna_link=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/gencode.vM38.transcripts.fa.gz
# gtf_link=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/gencode.vM38.basic.annotation.gtf.gz

## fly, chromosome: 2L, 2R, 3L, 3R, 4, X, Y, mitochondrion_genome, ...xxx
fasta_link=https://ftp.ensembl.org/pub/release-116/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.54.dna.toplevel.fa.gz
gtf_link=https://ftp.ensembl.org/pub/release-116/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.54.63.gtf.gz
cdna_link=https://ftp.ensembl.org/pub/release-116/fasta/drosophila_melanogaster/cdna/Drosophila_melanogaster.BDGP6.54.cdna.all.fa.gz

## cel, chromosome: I, II, III, IV, V, X, MtDNA
# fasta_link=https://ftp.ensembl.org/pub/release-116/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz
# gtf_link=https://ftp.ensembl.org/pub/release-116/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.63.gtf.gz
# cdna_link=https://ftp.ensembl.org/pub/release-116/fasta/caenorhabditis_elegans/cdna/Caenorhabditis_elegans.WBcel235.cdna.all.fa.gz


## set directory
# target_dir=/mnt/d/Index/hs/v49
target_dir=/mnt/d/Index/fly
root_dir=$target_dir && mkdir -p $root_dir && cd $root_dir && echo ">>>>>> start to processing <<<<<<<"1


#################################  process data  #################################
# fasta
fasta_gz=$(basename "$fasta_link")
fasta=${fasta_gz%.gz}
[ ! -f "$fasta_gz" ] && echo "downloading fasta" && wget $fasta_link -O $fasta_gz
[ ! -f "${fasta}.fai" ] && gunzip -k $fasta_gz && samtools faidx $fasta

# gtf
gtf_gz=$(basename "$gtf_link")
gtf=${gtf_gz%.gz}
[ ! -f "$gtf_gz" ] && echo "downloading gtf" && wget $gtf_link -O $gtf_gz
[ ! -f "$gtf" ] && gunzip -k $gtf_gz

# cdna
if [ -z "${cdna_link+x}" ]; then
    echo "cdna_link does not exist, extract cdna first"
    cdna_gz=cdna.fa.gz
	gffread -w - -g $fasta_gz $gtf_gz | gzip > "$cdna_gz"
else
    cdna_gz=$(basename "$cdna_link")
	[ ! -f "$cdna_gz" ] && echo "downloading cdna" && wget $cdna_link -O $cdna_gz
fi



#################################  salmon index  #################################
# salmon use ultra fast quasi mapping on transcript instead of genome
# mkdir -p salmon
# generate decoy
# zgrep "^>" $fasta_gz | cut -d " " -f 1 > decoys.txt && sed -i -e 's/>//g' decoys.txt
# add custome genome, usually like FPs, Cas9, resistance genes
# if [[ -f "custom_genome.fa" && -f "custom_genome.gtf" ]]; then
	# echo "add custome genome ...."
	# gffread -w custom_cdna.fa -g custom_genome.fa custom_genome.gtf && gzip custom_cdna.fa
	# gzip -k custom_genome.fa && \
		# zgrep "^>" custom_genome.fa.gz | cut -d " " -f 1 > custom_decoys.txt && \
		# sed -i -e 's/>//g' custom_decoys.txt && cat custom_decoys.txt >> decoys.txt
	# generate gentrome
	# zcat custom_cdna.fa.gz $cdna_gz custom_genome.fa.gz $fasta_gz > gentrome.fa.gz
# else
	# zcat $cdna_gz $fasta_gz > gentrome.fa.gz
# fi
# index
# salmon index -p 8 -t gentrome.fa.gz -d decoys.txt -i salmon --gencode
# remove unused files
# if [ $? -eq 0 ]; then rm custom_decoys.txt custom_genome.fa.gz custom_cdna.fa.gz gentrome.fa.gz custom_genome.fa.fai; fi


#################################  bowtie2 index  #################################
# note that bowtie2 is not a junction aware aligner
mkdir -p bowtie2
bowtie2-build --threads 8 $fasta_gz bowtie2/bowtie2


#################################  star index  #################################
# star is specifically designed for RNA mapping on genome
# mkdir -p star
# STAR --runThreadN 8 --runMode genomeGenerate --genomeDir star \
	# --genomeFastaFiles $fasta --sjdbGTFfile $gtf \
	# --sjdbOverhang 149 --limitGenomeGenerateRAM 30000000000 # ~ x/1024/1000/1000 GB
