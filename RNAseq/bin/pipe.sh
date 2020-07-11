#! /bin/sh

## read user_config.yml
source bin/yaml.sh
create_variables bin/config.yml


## comfirm info, if error, script will not excute, handle exits from shell or function but don't exit interactive shell
echo
echo "########## Please comfirm alignment information below ##########"
echo -e "Working directory: $work_dir\nSpecies: $species\nSequencing type: $seq_type\nFASTQ sufffix: $fastq_suffix\nThreads: $thread\nOutput_dir: output_dir"
echo "########## Please comfirm alignment information above ##########"
read -p "Choose Y/y to confirm: " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
  [[ "$0" = "$BASH_SOURCE" ]] && exit 1 || return 1 
fi


## load extdata file
if [ "$species" == "mm" ]; then
  star_index=$star_index_dir_mm
  salmon_index=$salmon_index_dir_mm
  gtf_file=$gtf_file_mm
  genome_size=2652783500
else
  star_index=$star_index_dir_hs
  salmon_index=$salmon_index_dir_hs
  gtf_file=$gtf_file_hs
  genome_size=2913022398
fi

## strandnesss for featurecounts
if [ "$lib_type" == "none" ]; then
  strandness=0
else
  if [ "$lib_type" == "reverse" ]; then
    strandness=2
  else
    strandness=1
  fi
fi




#------------------------------------------------------------------------------------------------------------------------

##### fastq fastqc #####
source bin/fastqc.sh $work_dir/$fastq_dir $work_dir/$fastqc_dir $fastq_suffix $thread
echo ""
sleep 1



##### trim fastq use fastp #####
source bin/fastp.sh $work_dir/$fastq_dir $work_dir/$trim_dir $fastq_suffix $read1_name $read2_name $thread
echo ""
sleep 1



##### salmon #####
source bin/salmon.sh $work_dir/$trim_dir $work_dir/$salmon_dir $read1_name $read2_name $thread $salmon_index
echo ""
sleep 1



##### star #####
if [ "$output_bam" == "true" ]; then
  source bin/star.sh $work_dir/$trim_dir $work_dir/$star_dir $read1_name $read2_name $thread $star_index
  echo ""
  sleep 1
fi

##### deeptools #####
if [ "$output_bam" == "true" ]; then
  source bin/deeptools.sh $work_dir/$star_dir $genome_size
  echo ""
  sleep 1
fi

##### feacturecounts #####
if [ "$output_bam" == "true" ]; then
  source bin/featurecounts.sh $work_dir/$star_dir $work_dir/$featurecounts_dir $thread $strandness $gtf_file
  echo ""
  sleep 1
fi



##### multiQC #####
source bin/multiqc.sh $work_dir $work_dir/$multiqc_dir
echo ""



##### calculate expression #####
#echo "============================"
#echo "gene expression quantification ..."
#Rscript bin/expr.R
#echo "============================"
#echo ""



#------------------------------------------------------------------------------------------------------------------------
echo "done!"











