#!/usr/bin/env bash
# MAGeCK CRISPR screen analysis (trimming -> count -> RRA -> MLE)
# used under Ubuntu 22.04, WSL2

set -euo pipefail

# check barcode manually
# zless xxx.fastq.gz

# conda activate mageck

project_dir="${1:-/media/hao/Data/Others/Tie2_SL/Seq_Raw_Data/w_rep}"
rawfastq_subdir="rawdata"
cleanfastq_subdir="clean_fastq"
counts_dir="counts"
rra_dir="RRA_analysis"
mle_dir="MLE_analysis"
library_file="${2:-/media/hao/Data1/CRISPR_screen/Library_GeCKOv2/metadata/Human_GeCKOv2_Library_combine.csv}"



cd "$project_dir" || { echo "ERROR: cannot cd to $project_dir" >&2; exit 1; }


echo "======================================================="
echo "================== cutadpat trimming =================="
echo "======================================================="
mkdir -p "$cleanfastq_subdir"

# for file in "$rawfastq_subdir"/**/*.fq.gz;
for file in "$rawfastq_subdir"/*.fq.gz;
do
	# only process read1 (for CRISPR use only)
	if [[ "$file" =~ _([rR])?1.f(ast)?q.gz ]]; then
	    fname=$(basename "$file")
	    echo -e "\n>>>>> cutadapt on $fname ..."
	    cutadapt -j 2 -Z -m 20 --discard-untrimmed \
	    	-g cttgtggaaaggacgaaacaccg \
	    	-o "$cleanfastq_subdir/$fname" \
	    	"$file"
	fi
done


echo "======================================================="
echo "==================== fastq prasing ===================="
echo "======================================================="
# get all gzipped fastq files
mapfile -t fastq_files < <(find "$cleanfastq_subdir" -type f -name "*.gz" | sort -u)
# echo "${fastq_files[@]}"
# printf '%s\n' "${fastq_files[@]}"

# get fastq file suffix
suffix_read1=$(printf '%s\n' "${fastq_files[@]}" | grep -o "_([rR])?1.f(ast)?q.gz" | sort -u)
echo -e "\nextract ${suffix_read1} as read-1 suffix"

# get fastq file read 1/2
fastq_read1=($(printf '%s\n' "${fastq_files[@]}" | grep "${suffix_read1}" | sort -u))

# get fastq file basenames
basenames=()
for i in "${fastq_read1[@]}"
	do
		name=$(basename "$i")
		basename=${name%${suffix_read1}}
		basenames+=("$basename")
	done
echo -e "\n------------- pharse samples ------------"
printf '%s\n' "${basenames[@]}"
echo -e "-------------------------------------------------\n"



echo "======================================================="
echo "===================== MaGeCK count ===================="
echo "======================================================="
mkdir -p "$counts_dir"

# get full paths
full_paths=("${basenames[@]/#/$cleanfastq_subdir/}")
full_paths=("${full_paths[@]/%/$suffix_read1}")
# printf '%s\n' "${full_paths[@]}"

# mageck count -l $library_file -n $counts_dir/mageck_count --pdf-report --test-run\
#	--sample-label "s1,s2" \
#	--fastq clean_fastq/R10_CKDL240002461-1A_H2NMCDSXC_L2_1.fq.gz clean_fastq/R2_CKDL240002453-1A_H2NMCDSXC_L2_1.fq.gz

mageck count -l "$library_file" -n "$counts_dir/mageck" \
	--sample-label "$(IFS=,; printf '%s' "${basenames[*]}")" \
	--fastq "${full_paths[@]}"



echo "======================================================="
echo "===================== MaGeCK RRA ======================"
echo "======================================================="
mkdir -p "$rra_dir"

compare_name=(
	"shNT"
	"shTEK"
	"shTEK_shNT"
	)
compare_t=(
	"N71_H32Y5ALXX_L1_rep1,N71_H32Y5ALXX_L1_rep2"
	"N73_H32Y5ALXX_L1_rep1,N73_H32Y5ALXX_L1_rep2"
	"N73_H32Y5ALXX_L1_rep1,N73_H32Y5ALXX_L1_rep2"
	)
compare_c=(
	"NC_H32Y5ALXX_L1_rep1,NC_H32Y5ALXX_L1_rep2"
	"NC_H32Y5ALXX_L1_rep1,NC_H32Y5ALXX_L1_rep2"
	"N71_H32Y5ALXX_L1_rep1,N71_H32Y5ALXX_L1_rep2"
	)

for i in "${!compare_name[@]}"; do
  # printf "%s\t%s\n" "$i" "${compare_name[$i]}"
  printf "%s\t%s\n" "$i" "${compare_name[$i]}"
  mageck test -k "$counts_dir/mageck.count.txt" \
  	-n "$rra_dir/${compare_name[$i]}" \
		-t "${compare_t[$i]}" \
		-c "${compare_c[$i]}"
done



echo "======================================================="
echo "===================== MaGeCK MLE ======================"
echo "======================================================="
mkdir -p "$mle_dir"

mageck mle -k "$counts_dir/mageck.count.txt" \
 -d MAGeCK_mle_design.txt -n "$mle_dir/mle"


# mageck mle -k $counts_dir/mageck.count.txt\
#  -d MAGeCK_mle_design_APP.txt -n $mle_dir/APP

#  mageck mle -k $counts_dir/mageck.count.txt\
#  -d MAGeCK_mle_design_NOTCH1.txt -n $mle_dir/NOTCH1
 
#  mageck mle -k $counts_dir/mageck.count.txt\
#  -d MAGeCK_mle_design_Gal4DBD.txt -n $mle_dir/Gal4DBD
