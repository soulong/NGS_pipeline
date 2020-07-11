#! /bin/sh


input_dir=$1
output_dir=$2
fastq_suffix=$3
read1_name=$4
read2_name=$5
thread=$6



mkdir -p $output_dir 
echo "============================"

echo "quality trim by fastp ..."
echo `date` >> $output_dir/fastp.log

filelist_uni=$(find $input_dir -type f -name "*${read1_name}${fastq_suffix}")



for i in $filelist_uni
    do
      
      fullname=$(basename $i)
      fullpath=$(dirname $i)
      basename=${fullname%${read1_name}${fastq_suffix}}


      if [ -f "$output_dir/$basename${read1_name}.fastq.gz" ]; then
        echo "$basename exists, jump to next one"
        continue
      fi

      trap 'echo stopped, removing created files; \
        rm $output_dir/$basename*; \
        exit 1' SIGINT SIGTERM

      echo "trim on $basename ..."

      SECONDS=0

      fastp -i ${fullpath}/${basename}${read1_name}${fastq_suffix} \
          -I ${fullpath}/${basename}${read2_name}${fastq_suffix} \
          -o ${output_dir}/${basename}${read1_name}.fastq.gz \
          -O ${output_dir}/${basename}${read2_name}.fastq.gz \
          -l 36 -j ${output_dir}/${basename}.json -h ${output_dir}/${basename}.html \
          -R "$basename" -w $thread \
          2>> ${output_dir}/fastp.log
          
      duration=$SECONDS
      echo ">> $(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed"
      sleep 1
  done

trap '' SIGINT SIGTERM

sleep 1

echo "============================"
