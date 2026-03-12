#!/usr/bin/env bash
# helper.sh - Common functions for genomic pipelines
# Usage: source this file in your main pipeline script

set -euo pipefail

# ============================= Logging Functions =============================
# Usage: log "message"
log() { echo "[$(date '+%Y-%m-%d %H:%M')] $*" >&2; }
log_info() { log "[INFO] $*"; }
log_warn() { log "[WARN] $*"; }
log_error() { log "[ERROR] $*"; }

# Usage: die "error message"
die() { log_error "$*"; exit 1; }


# ============================= File Validation =============================
require_dirs() {
    # Usage: require_dirs <path1> [<path2> ...] or "path1 path2 ..."
    local -a paths=()
    for a in "$@"; do
        # split whitespace-separated strings into individual tokens
        read -ra parts <<< "$a"
        for p in "${parts[@]}"; do
            paths+=("$p")
        done
    done

    for p in "${paths[@]}"; do
        if [[ ! -d "$p" ]]; then
            die "Missing required directory: $p"
        fi
    done
}

require_files() {
    # Usage: require_files <file1> [<file2> ...] or "file1 file2 ..."
    local -a files=()
    for a in "$@"; do
        # split whitespace-separated strings into individual tokens
        read -ra parts <<< "$a"
        for p in "${parts[@]}"; do
            files+=("$p")
        done
    done

    for f in "${files[@]}"; do
        if [[ ! -f "$f" ]]; then
            die "Missing required file: $f"
        fi
    done
}

require_commands() {
    # Usage: require_commands <cmd1> [<cmd2> ...] or "cmd1 cmd2 ..."
    local -a cmds=()
    for a in "$@"; do
        # split whitespace-separated strings into individual tokens
        read -ra parts <<< "$a"
        for p in "${parts[@]}"; do
            cmds+=("$p")
        done
    done

    for c in "${cmds[@]}"; do
        # only check non-empty tokens
        [[ -z "$c" ]] && continue
        if ! command -v "$c" &>/dev/null; then
            die "Required command not found: $c"
        fi
    done
}

create_dirs() {
    # Usage: safe_create_dirs <directory1> [<directory2> ...]
    for dir in "$@"; do
        if [[ ! -d "$dir" ]]; then
            if mkdir -p "$dir"; then
                log_info "Created directory: $dir"
            else
                die "Failed to create directory: $dir"
            fi
        fi
    done
}

check_files_exists() {
    # Usage: check_files_exists <path1> [<path2> ...]
    # Returns: 0 if ALL exist and valid, 1 otherwise
    local paths=("$@")
    
    for path in "${paths[@]}"; do
        if [[ ! -e "$path" ]]; then
            return 1
        fi
        
        if [[ -f "$path" ]]; then
            local size
            size=$(stat -f%z "$path" 2>/dev/null || stat -c%s "$path" 2>/dev/null || echo "0")
            if [[ $size -lt 1 ]]; then
                return 1
            fi
        fi
    done
    
    return 0
}

get_base_fq_name() {
    # Extract base name from FASTQ file, removing common extensions
    # Usage: get_base_name "sample_R1.fastq.gz"
    # Returns: sample_R1
    local file="$1"
    local name
    name=$(basename "$file")
    name=${name%.fastq.gz}
    name=${name%.fq.gz}
    name=${name%.fastq}
    name=${name%.fq}
    echo "$name"
}

# ============================= Process Management =============================
setup_trap() {
    # Setup trap for interrupt signals with cleanup
    # Usage: setup_trap <cleanup_pattern>
    # Example: setup_trap "output_dir/sample*"
    local cleanup_pattern="$1"
    
    trap "echo 'Interrupted! Cleaning up...'; rm -rf $cleanup_pattern; exit 1" SIGINT SIGTERM
}

clear_trap() {
    # Clear interrupt traps, combinated with SIGINT and SIGTERM (with setup_trap more often)
    # Usage: clear_trap
    trap '' SIGINT SIGTERM
}

# ============================= File Processing =============================
loop_over_csv() {
    # Safely read CSV file, handling DOS line endings and headers
    # Usage: loop_over_csv <csv_file> <callback_function>
    # Callback receives: sample group control target fq1 fq2
    local file="$1"
    local callback="$2"
    
    while IFS=, read -r sample group control target fq1 fq2 || [[ -n "${sample:-}" ]]; do
        # Skip header
        [[ "$sample" == "sample" ]] && continue
        # Skip empty lines
        [[ -z "$sample" ]] && continue
        
        "$callback" "$sample" "$group" "$control" "$target" "$fq1" "$fq2"
    done < <(tail -n +2 "$file" | tr -d '\r')
}

generate_samplesheet() {
    # Auto-generate samplesheet by scanning for paired-end FASTQ files
    # Usage: generate_samplesheet [input_dir] [output_csv]
    # Assumes R1/R2 naming convention
    local input_dir="${1:-.}"
    local output_csv="${2:-samplesheet.csv}"
    
    if [[ ! -d "$input_dir" ]]; then
        log_error "Directory does not exist: $input_dir"
        return 1
    fi
    
    log_info "Generating samplesheet from: $input_dir"
    
    # Write header
    printf "sample,group,control,target,fq1,fq2\n" > "$output_csv"
    
    # Find R1 files and infer R2
    find "$input_dir" -type f -regextype posix-extended \
        -regex '.*_[rR]?1(_.*)?\.f(ast)?q\.gz' | sort | while IFS= read -r fq1; do
        
        local fq1_dir fq1_base fq2_base fq2
        fq1_dir=$(dirname "$fq1")
        fq1_base=$(basename "$fq1")
        
        # Convert R1 to R2
        fq2_base=$(echo "$fq1_base" | sed -E 's/^(.*)(_[rR]?)1((_.*)?\.f(ast)?q\.gz)$/\1\22\3/')
        fq2="$fq1_dir/$fq2_base"
        
        if [[ ! -f "$fq2" ]]; then
            log_warn "Missing mate pair: $fq2 (derived from $fq1)"
            fq2="NA"
        fi
        
        # Extract sample name by removing common suffixes
        local sample_name
        sample_name=$(basename "$fq1" | sed -E 's/_R?1(_.*)?\.f(ast)?q\.gz$//')
        
        printf "%s,,,,%s,%s\n" "$sample_name" "$fq1" "$fq2" >> "$output_csv"
    done
    
    log_info "Generated samplesheet: $output_csv"
}



# ============================= Statistics Extraction =============================
extract_flagstat_value() {
    # Extract value from samtools flagstat output
    # Usage: extract_flagstat_value <flagstat_file> <pattern> [column]
    # Example: extract_flagstat_value "sample.flagstat" "mapped" 1
    local file="$1"
    local pattern="$2"
    local col="${3:-1}"
    
    if [[ ! -f "$file" ]]; then
        log_warn "Flagstat file not found: $file"
        echo "0"
        return 1
    fi
    
    awk -v p="$pattern" -v c="$col" '
        $0 ~ p { print $c; exit }
    ' "$file"
}

extract_salmon_stats() {
    # Extract mapping statistics from Salmon log
    # Usage: extract_salmon_stats <salmon_dir/sample>
    # Returns: "total_reads mapped_reads mapping_rate"
    local salmon_dir="$1"
    local log_file="$salmon_dir/logs/salmon_quant.log"
    
    if [[ ! -f "$log_file" ]]; then
        echo "0 0 0.0"
        return 1
    fi
    
    local total_reads mapped_reads mapping_rate
    
    # Extract total fragments processed
    total_reads=$(grep -oP "Observed \K[0-9,]+" "$log_file" | tr -d ',' || echo "0")
    
    # Extract mapping rate
    mapping_rate=$(grep -oP "Mapping rate = \K[0-9.]+%" "$log_file" | tr -d '%' || echo "0.0")
    
    # Calculate mapped reads
    mapped_reads=$(awk -v total="$total_reads" -v rate="$mapping_rate" \
        'BEGIN{printf "%.0f", total * rate / 100}')
    
    echo "$total_reads $mapped_reads $mapping_rate"
}

extract_star_stats() {
    # Extract mapping statistics from STAR log
    # Usage: extract_star_stats <star_log_file>
    # Returns: "total_reads uniquely_mapped multi_mapped unmapped"
    local log_file="$1"
    
    if [[ ! -f "$log_file" ]]; then
        echo "0 0 0 0"
        return 1
    fi
    
    local total uniquely multi unmapped
    
    total=$(grep "Number of input reads" "$log_file" | awk '{print $NF}')
    uniquely=$(grep "Uniquely mapped reads number" "$log_file" | awk '{print $NF}')
    multi=$(grep "Number of reads mapped to multiple loci" "$log_file" | awk '{print $NF}')
    unmapped=$(grep "Number of reads unmapped: too short" "$log_file" | awk '{print $NF}')
    
    echo "${total:-0} ${uniquely:-0} ${multi:-0} ${unmapped:-0}"
}

# ============================= Common Analysis Functions =============================
# Pre-process sample, run fastqc and fastp trimming
# generally loop over samples
run_preprocess() {
    # sample="$1" group="$2" control="$3" target="$4" fq1="$5" fq2="$6"
    local sample="$1"
    local r1="$5"
    local r2="$6"
    
    # # FastQC
    # local qc1="$FASTQC_DIR/${sample}_R1_fastqc.zip"
    # local qc2="$FASTQC_DIR/${sample}_R2_fastqc.zip"
    # if check_files_exists "$qc1" "$qc2"; then
    #     log_info "[$sample] FastQC completed"
    # else
    #     log_info "[$sample] Running FastQC"
    #     fastqc -t "$THREADS" -o "$FASTQC_DIR" --quiet "$r1" "$r2" 2>/dev/null || {
    #         log_error "[$sample] FastQC failed"
    #         return 1
    #     }
    #     mv "$FASTQC_DIR/$(get_base_fq_name $r1)_fastqc.html" "$FASTQC_DIR/${sample}_R1_fastqc.html"
    #     mv "$FASTQC_DIR/$(get_base_fq_name $r1)_fastqc.zip" "$qc1"
    #     mv "$FASTQC_DIR/$(get_base_fq_name $r2)_fastqc.html" "$FASTQC_DIR/${sample}_R2_fastqc.html"
    #     mv "$FASTQC_DIR/$(get_base_fq_name $r2)_fastqc.zip" "$qc2"
    # fi

    # fastp trimming
    local clean_r1="$FASTP_DIR/${sample}_R1.fq.gz"
    local clean_r2="$FASTP_DIR/${sample}_R2.fq.gz"
    if check_files_exists "$clean_r1" "$clean_r2"; then
        log_info "[$sample] Trimming completed"
    else
        log_info "[$sample] Running fastp trim"
        setup_trap "$FASTP_DIR/${sample}*"
        fastp -w "$THREADS" -i "$r1" -I "$r2" -o "$clean_r1" -O "$clean_r2" \
            -j "$FASTP_DIR/${sample}.json" -h "$FASTP_DIR/${sample}.html" \
            -R "$sample" 2>> "$FASTP_DIR/fastp.log" || { 
            log_error "[$sample] fastp failed"
            clear_trap
            return 1
        }
        clear_trap
    fi
}


# ============================= Normalization =============================
calculate_scale_factor() {
    # Calculate normalization scale factor
    # Usage: calculate_scale_factor <method> <filtered_reads> [spike_reads] [sample] [spike_free_file]
    # Returns: scale factor as float
    local method="$1"
    local filtered_reads="$2"
    local spike_reads="${3:-0}"
    local sample="${4:-unknown}"
    local spike_free_file="${5:-}"

    local scale_factor="1.0"
    
    case "$method" in
        CPM)
            if [[ $filtered_reads -gt 0 ]]; then
                scale_factor=$(awk "BEGIN{printf \"%.6f\", 1e6 / $filtered_reads}")
            else
                log_warn "Filtered reads = 0, using scale factor 1.0"
            fi
            ;;
        Spike)
            if [[ $spike_reads -gt 0 ]]; then
                scale_factor=$(awk "BEGIN{printf \"%.6f\", 1e4 / $spike_reads}")
            else
                log_warn "Spike reads = 0, using scale factor 1.0"
            fi
            ;;
        SpikeFree)
            if [[ -z "$spike_free_file" || ! -f "$spike_free_file" ]]; then
                die "SpikeFree method requires spike-free scale factor file"
            fi
            
            local sf_raw
            sf_raw=$(awk -F'\t' -v s="$sample" '
                NR>1 {
                    gsub(/\.filtered\.bam$/, "", $1)
                    if ($1 == s) { print $7; exit }
                }
            ' "$spike_free_file")
            
            if [[ -z "$sf_raw" ]]; then
                die "Sample $sample not found in $spike_free_file"
            fi
            
            scale_factor=$(awk -v fr="$filtered_reads" -v sf="$sf_raw" \
                'BEGIN{printf "%.6f", 1e6 / (fr * sf)}')
            ;;
        *)
            log_error "Unknown normalization method: $method"
            return 1
            ;;
    esac
    
    echo "$scale_factor"
}

# ============================= Peak Processing =============================
merge_peakfiles() {
    # Merge multiple peak files into consensus BED
    # Usage: merge_peakfiles --input-peaks file1 file2 --output-bed output.bed --chrom-size chrom.sizes
    local input_peaks=""
    local output_bed=""
    local chrom_size=""

    while [[ $# -gt 0 ]]; do
        case "$1" in
            --input-peaks)
                input_peaks="$2"
                shift 2
                ;;
            --output-bed)
                output_bed="$2"
                shift 2
                ;;
            --chrom-size)
                chrom_size="$2"
                shift 2
                ;;
            -h|--help)
                echo "Usage: make_consensus_peak --input-peaks file1 file2 --output-bed output.bed --chrom-size chrom.sizes"
                echo "  --input-peaks    Space-separated list of peak files"
                echo "  --output-bed     Output consensus BED file"
                echo "  --chrom-size     Chromosome sizes file for sorting"
                return 0
                ;;
            *)
                log_error "Unknown option: $1"
                return 1
                ;;
        esac
    done

    # Validate inputs
    [[ -z "$input_peaks" ]] && { log_error "Input peaks required (--input-peaks)"; return 1; }
    [[ -z "$output_bed" ]] && { log_error "Output file required (--output-bed)"; return 1; }
    [[ -z "$chrom_size" ]] && { log_error "Chromosome sizes file required (--chrom-size)"; return 1; }
    require_files "$chrom_size"

    # Convert space-separated string to array
    local -a peak_files
    read -ra peak_files <<< "$input_peaks"

    if [[ ${#peak_files[@]} -eq 0 ]]; then
        log_error "No peak files provided"
        return 1
    fi

    # Validate all peak files exist
    local missing=0
    for pf in "${peak_files[@]}"; do
        if [[ ! -f "$pf" ]]; then
            log_warn "Peak file not found: $pf"
            ((missing++))
        fi
    done

    if [[ $missing -gt 0 ]]; then
        log_error "Some peak files are missing"
        return 1
    fi

    # log_info "Merging ${#peak_files[@]} peak files into $output_bed"

    # Merge peaks: concatenate -> extract coords -> sort -> merge -> sort by genome
    cat "${peak_files[@]}" | \
        cut -f1-3 | \
        bedtools sort 2>/dev/null | \
        bedtools merge 2>/dev/null | \
        bedtools sort -g "$chrom_size" 2>/dev/null > "$output_bed"

    if [[ -s "$output_bed" ]]; then
        local count
        count=$(wc -l < "$output_bed")
        log_info "$count regions merged"
        return 0
    else
        log_error "Failed to generate consensus peaks"
        rm -f "$output_bed"
        return 1
    fi
}


# ============================= Genomic Indices Setup =============================
setup_genome() {
    # Setup genome indices for either CutTag or RNAseq based on species
    # Usage: setup_genome <species> <index_rootdir>
    local species="$1"
    local index_rootdir="$2"
    
    log_info "Setting up genome indices for species: $species"
 
    case "$species" in
        chm13)
            INDEX_DIR="$index_rootdir/hs/chm13"
            FASTA="${INDEX_DIR}/chm13v2.0.fa.gz"
            GTF="${INDEX_DIR}/chm13v2.0.gtf.gz"
            CHROM_SIZES="${FASTA}.fai"
            GSIZE=2913022398
            ;;
        hs)
            INDEX_DIR="$index_rootdir/hs/v49"
            FASTA="${INDEX_DIR}/GRCh38.primary_assembly.genome.fa.gz"
            GTF="${INDEX_DIR}/gencode.v49.basic.annotation.gtf.gz"
            CHROM_SIZES="${FASTA}.fai"
            GSIZE=2913022398
            ;;
        mm)
            INDEX_DIR="$index_rootdir/mm/vM38"
            FASTA="${INDEX_DIR}/GRCm39.primary_assembly.genome.fa.gz"
            GTF="${INDEX_DIR}/gencode.vM38.annotation.gtf.gz"
            CHROM_SIZES="${FASTA}.fai"
            GSIZE=2654621783
            ;;
    esac
    
    # Setup common indices
    SPIKE_INDEX="$index_rootdir/Ecoli_novoprotein/bowtie2/bowtie2"

    BOWTIE2_INDEX="$INDEX_DIR/bowtie2/bowtie2"
    SALMON_INDEX="$INDEX_DIR/salmon"
    STAR_INDEX="$INDEX_DIR/star"
        
    # Setup default BED regions for peak heatmaps
    BED_REGIONS=(
        "$INDEX_DIR/genes.bed"
        "$INDEX_DIR/genes_protein_coding.bed"
    )
}



# ============================= Merge fq by same samples =============================
# Merges split FASTQ files (from SRA / multiple lanes/runs) into one pair per sample
# Used for downloaded from GEO/ENA, example:
# mamba activate ngs && fastq-dl --group-by-sample --ignore --outdir fastq --cpus 16 --force -a SRP108500
# Arguments:
#   $1 = path to samplesheet (SraRunTable.csv, downloaded from SRA selector metadata)
#   $2 = directory containing the raw *_1.fastq.gz and *_2.fastq.gz files (fastq-dl style)
#   $3 = output directory (will be created if missing)   [default: ./merged]
# Assumptions (this is SRA selector metadata style):
#   - Column 1  = Run (SRRxxxxxx)
#   - Column 25 = Sample Name (usually GSMxxxxxx)
#   - Files are named exactly like: SRR5638491_1.fastq.gz / SRR5638491_2.fastq.gz
# =============================================================================
merge_fastq_by_sample() {
    local samplesheet="$1"
    local input_dir="$2"
    local out_dir="${3:-merged}"

    if [[ ! -f "$samplesheet" ]]; then
        echo "Error: samplesheet not found: $samplesheet" >&2
        return 1
    fi
    if [[ ! -d "$input_dir" ]]; then
        echo "Error: input directory not found: $input_dir" >&2
        return 1
    fi

    mkdir -p "$out_dir" || return 1

    # Extract unique sample names (skip header)
    local samples
    samples=$(awk -F',' 'NR>1 {print $25}' "$samplesheet" | sort -u)

    local total=$(echo "$samples" | wc -l)
    local i=0

    echo "Found $total unique samples in $samplesheet"

    while IFS= read -r sample; do
        ((i++))
        printf "\n[%2d/%d] %s\n" "$i" "$total" "$sample"

        # Find all runs for this sample
        local runs
        runs=$(awk -F',' -v s="$sample" '$25==s {print $1}' "$samplesheet")

        # Build list of R1 and R2 files (only if they exist)
        local r1_files=()
        local r2_files=()

        for run in $runs; do
            local f1="${input_dir}/${run}_1.fastq.gz"
            local f2="${input_dir}/${run}_2.fastq.gz"

            [[ -f "$f1" ]] && r1_files+=("$f1")
            [[ -f "$f2" ]] && r2_files+=("$f2")
        done

        if [[ ${#r1_files[@]} -eq 0 ]]; then
            echo "  Warning: no R1 files found for $sample" >&2
            continue
        fi

        # Merge R1
        echo "  Merging ${#r1_files[@]} R1 files → ${out_dir}/${sample}_1.fastq.gz"
        cat "${r1_files[@]}" > "${out_dir}/${sample}_1.fastq.gz"

        # Merge R2 (if any exist)
        if [[ ${#r2_files[@]} -gt 0 ]]; then
            if [[ ${#r1_files[@]} -ne ${#r2_files[@]} ]]; then
                echo "  Warning: number of R1 (${#r1_files[@]}) ≠ R2 (${#r2_files[@]}) for $sample" >&2
            fi
            echo "  Merging ${#r2_files[@]} R2 files → ${out_dir}/${sample}_2.fastq.gz"
            cat "${r2_files[@]}" > "${out_dir}/${sample}_2.fastq.gz"
        else
            echo "  Note: no R2 files found for $sample" >&2
        fi

    done <<< "$samples"

    echo -e "\nDone. Merged files are in: $out_dir/"
}