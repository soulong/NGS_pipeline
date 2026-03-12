#! /bin/bash
# ubuntu_24.04_LTS
# last edited at 2025-11-13 by Hao He

#################################  setup  #################################
conda activate ngs
# conda env export --no-builds | grep -v "^prefix: " > ngs.yml
# conda env create -f ngs.yml

# set directory
root_dir=/mnt/d/Index/hs/chm13 \
	&& mkdir -p $root_dir && cd $root_dir \
	&& echo ">>>>>> start to processing <<<<<<<"


## link can be remote download link or local file
# t2t-chm13
fasta_link="chm13v2.0.fa.gz"
gtf_link="chm13v2.0.gtf.gz"

# # hs
# fasta_link=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.primary_assembly.genome.fa.gz
# gtf_link=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.basic.annotation.gtf.gz

# # mm
# fasta_link=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/GRCm39.primary_assembly.genome.fa.gz
# gtf_link=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/gencode.vM38.basic.annotation.gtf.gz

# get filename from link
fasta_gz=$(basename "$fasta_link")
gtf_gz=$(basename "$gtf_link")


#################################  download and process data  #################################
# bgzip -@ 4 -k $fasta_gz
# file $fasta_gz

## fasta
[ ! -f "$fasta_gz" ] && echo "downloading fasta" && wget $fasta_link -O $fasta_gz
zcat "$fasta_gz" | bgzip -@ 4 -c > "${fasta_gz}.tmp" && mv "${fasta_gz}.tmp" "$fasta_gz"
[ ! -f "${fasta_gz}.fai" ] && echo "index fasta" && samtools faidx -@ 4 $fasta_gz

## gtf
[ ! -f "$gtf_gz" ] && echo "downloading gtf" && wget $gtf_link -O $gtf_gz
zcat "$gtf_gz" | bgzip -@ 4 -c > "${gtf_gz}.tmp" && mv "${gtf_gz}.tmp" "$gtf_gz"

# get ungzipped files
[[ ! -s "${fasta_gz%.gz}" ]] && gzip -d -k $fasta_gz && samtools faidx -@ 4 ${fasta_gz%.gz}
[[ ! -s "${gtf_gz%.gz}" ]] && gzip -d -k $gtf_gz


#################################  salmon index  #################################
# salmon use ultra fast quasi mapping on transcript instead of genome
mkdir -p salmon
# generate decoy
zgrep "^>" $fasta_gz | cut -d " " -f 1 > decoys.txt && sed -i -e 's/>//g' decoys.txt
# generate cdna, generate don't support .gz files
zcat "$gtf_gz" | gffread - -g ${fasta_gz%.gz} -w cdna.fa

# add custome genome, usually like FPs, Cas9, resistance genes
custome_fa="../../chrN.fa"
custome_gtf="../../chrN.gtf"
if [[ -f "$custome_fa" && -f "$custome_gtf" ]]; then
	echo "add custome genome ...."
	gffread -g $custome_fa $custome_gtf -w custom_cdna.fa
	grep "^>" $custome_fa | cut -d " " -f 1 > custom_decoys.txt && \
		sed -i -e 's/>//g' custom_decoys.txt && \
		cat custom_decoys.txt >> decoys.txt
	# generate gentrome
	cat custom_cdna.fa cdna.fa $custome_fa ${fasta_gz%.gz} | \
		bgzip -@ 4 -c > gentrome.fa.gz
else
	cat cdna.fa ${fasta_gz%.gz} | bgzip -@ 4 -c > gentrome.fa.gz
fi

# index
salmon index -p 8 -t gentrome.fa.gz -d decoys.txt -i salmon --gencode
# remove unused files
if [ $? -eq 0 ]; then rm ${custome_fa}.fai custom_decoys.txt decoys.txt custom_cdna.fa cdna.fa gentrome.fa.gz; fi



#################################  bowtie2 index  #################################
# note that bowtie2 is not a junction aware aligner
mkdir -p bowtie2
bowtie2-build --threads 8 $fasta_gz bowtie2/bowtie2


#################################  star index  #################################
# star is specifically designed for RNA mapping on genome
# star don't support compressed .fa and .gtf files
mkdir -p star
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir star \
	--genomeFastaFiles "${fasta_gz%.gz}" --sjdbGTFfile "${gtf_gz%.gz}" \
	--sjdbOverhang 149 --limitGenomeGenerateRAM 24000000000 # ~ x/1024/1000/1000 GB
# remove unused files
if [ $? -eq 0 ]; then rm ${fasta_gz%.gz} ${fasta_gz%.gz}.fai ${gtf_gz%.gz}; fi