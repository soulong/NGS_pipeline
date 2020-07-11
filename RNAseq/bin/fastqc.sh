#! /bin/sh


input_dir=$1
output_dir=$2
fastq_suffix=$3
thread=$4



echo "============================"


if [ ! -d "$output_dir" ]; then

	mkdir -p $output_dir


	echo "start fastqc ..."

	
	echo `date` >> $output_dir/fastqc.log

	filelist=$(find $input_dir -type f -name "*$fastq_suffix")

	SECONDS=0

	fastqc -t $thread -o $output_dir $filelist \
		2>> $output_dir/fastqc.log

	duration=$SECONDS
	echo ">> $(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed"
else
  echo "existing fastqc directory found, skip fastqc step"
fi

echo "============================"