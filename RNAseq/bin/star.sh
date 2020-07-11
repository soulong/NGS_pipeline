#! /bin/sh


input_dir=$1
output_dir=$2
read1_name=$3
read2_name=$4
thread=$5
star_index=$6


mkdir -p $output_dir

echo "============================"

echo "star algin ..."
echo `date` >> $output_dir/star.log


for i in $input_dir/*$read1_name.fastq.gz
do
	fullname=`basename $i`
	basename=${fullname%$read1_name.fastq.gz}

	if [ -d "$output_dir/$basename" ]; then
        echo "$basename exists, jump to next one"
        continue
    fi

    trap 'echo stopped, removing created files; \
        rm -r $output_dir/$basename*; \
        exit 1' SIGINT SIGTERM

	mkdir -p $output_dir/$basename

	echo "star align on $basename ..."

	SECONDS=0

	STAR --runThreadN $thread --genomeDir $star_index --readFilesCommand zcat \
	--readFilesIn $input_dir/$basename$read1_name.fastq.gz $input_dir/$basename$read2_name.fastq.gz \
	--outFileNamePrefix $output_dir/$basename/ \
	--outSAMtype SAM

	sleep 1

	echo "samtools sort ..."
	samtools view -@ $thread -Sbh $output_dir/$basename/Aligned.out.sam | \
	samtools sort -@ $thread -m 2G -o $output_dir/$basename/$basename.bam

	echo "samtools index ..."
	samtools index $output_dir/$basename/$basename.bam


	rm $output_dir/$basename/*.sam
	
	duration=$SECONDS
    echo ">> $(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed"
    sleep 1
    
done

trap '' SIGINT SIGTERM

sleep 1

echo "============================"






# echo "============================"

# echo "bedtools genomecov ..."


# bam_list=$(find $output_dir -type f -name "*.bam")

# for i in $bam_list
#     do
      
#       fullname=$(basename $i)
#       fullpath=$(dirname $i)
#       basename=${fullname%.bam}

#       if [ -f "$fullpath/${basename}.bedgraph.gz" ]; then
#         echo "${basename}.bedgraph.gz exists, jump to next one"
#         continue
#       fi

#       trap 'echo stopped, removing created files; \
#         rm $output_dir/$basename*; \
#         exit 1' SIGINT SIGTERM

#       echo "genomecov on $basename ..."

#       SECONDS=0

#       bedtools genomecov -ibam $i -bg \
#       | gzip > $fullpath/${basename}.bedgraph.gz

#       duration=$SECONDS
#       echo ">> $(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed"
#       echo ""
#       sleep 1
#   done

# trap '' SIGINT SIGTERM

# sleep 1


# echo "============================"





