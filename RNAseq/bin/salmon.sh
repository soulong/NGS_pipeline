#! /bin/sh


input_dir=$1
output_dir=$2
read1_name=$3
read2_name=$4
thread=$5
salmon_index=$6


mkdir -p $output_dir

echo "============================"

echo "salmon algin ..."
echo `date` >> $output_dir/salmon.log


for i in $input_dir/*${read1_name}.fastq.gz
do
	fullname=`basename $i`
	basename=${fullname%${read1_name}.fastq.gz}

	if [ -d "$output_dir/$basename" ]; then
        echo "$basename exists, jump to next one"
        continue
    fi

    trap 'echo stopped, removing created files; \
        rm -r $output_dir/$basename*; \
        exit 1' SIGINT SIGTERM
        
	echo ""
	echo "salmon quant on $basename ..."

	SECONDS=0

	salmon quant -i $salmon_index -l A \
        -1 $input_dir/$basename${read1_name}.fastq.gz -2 $input_dir/$basename${read2_name}.fastq.gz \
        -p $thread --seqBias --gcBias --validateMappings -o $output_dir/$basename \
		2>> $output_dir/salmon.log

	sleep 1

	duration=$SECONDS
    echo ">> $(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed"
    sleep 1
    
done

trap '' SIGINT SIGTERM

sleep 1

echo "============================"
