#! /bin/sh


input_dir=$1
output_dir=$2
thread=$3
strandness=$4
gtf_file=$5

mkdir -p $output_dir

echo "============================"
echo "featurecounts gene level quantification ..."

echo `date` >> $output_dir/featurecounts.log

if [ ! -f $output_dir/gene_featurecounts.txt ]; then
	
	SECONDS=0

	featureCounts -T $thread -s $strandness -t exon -g gene_id -Q 20 -p -C --donotsort -a $gtf_file \
		-o $output_dir/gene_featurecounts.txt `ls $input_dir/*/*.bam` \
		2>> $output_dir/featurecounts.log

	duration=$SECONDS
	echo ">> $(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed"
	sleep 1
else
	echo "featureCounts result has found, skip"
fi

echo "============================"
