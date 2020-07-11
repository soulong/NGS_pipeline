#! /bin/sh


bam_dir=$1
genome_size=$2



mkdir -p $output_dir


echo "============================"

echo "deeptools bedgraph and bigwig ..."


bam_list=$(find $bam_dir -type f -name "*.bam")

for i in $bam_list
    do
      
      fullname=$(basename $i)
      fullpath=$(dirname $i)
      basename=${fullname%.bam}

    trap 'echo stopped, removing created files; \
        rm $bam_dir/$basename/${basename}.bedgraph.gz $bam_dir/$basename/${basename}.bigwig; \
        exit 1' SIGINT SIGTERM


    SECONDS=0

	# if [ -f "$fullpath/${basename}.bedgraph" ]; then
 #        echo "${basename}.bedgraph.gz exists, jump to next one"

 #    else
	# 	echo "bedgraph on $basename ..."
 #      	bamCoverage --bam $i --normalizeTo1x $genome_size --binSize 10 -p max \
 #      		--outFileFormat bedgraph --minMappingQuality 20 --ignoreDuplicates \
 #      		-o $fullpath/${basename}.normlize1x.bedgraph
 #      	echo ""
 #    fi


	if [ -f "$fullpath/${basename}.normlize1x.bigwig" ]; then
        echo "${basename}.bigwig exists, jump to next one"

    else
	  echo "bigwig on $basename ..."
      bamCoverage --bam $i --normalizeTo1x $genome_size --binSize 10 -p max \
      	--outFileFormat bigwig --minMappingQuality 20 --ignoreDuplicates \
      	-o $fullpath/${basename}.normlize1x.bigwig
      	echo ""
    fi


    duration=$SECONDS
    echo ">> $(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed"
    echo ""

    sleep 1

  done

trap '' SIGINT SIGTERM

sleep 1


echo "============================"

