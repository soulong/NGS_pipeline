#! /bin/sh


input_dir=$1
output_dir=$2


echo "============================"
echo "multiqc report ..."

if [ ! -d $output_dir ]; then

	mkdir -p $output_dir

	SECONDS=0

	multiqc -f $input_dir -o $output_dir

	duration=$SECONDS
	echo ">> $(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed"
	sleep 1

else
	echo "multiqc directory has found, skip"
	sleep 1
fi

echo "============================"
