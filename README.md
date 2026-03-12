# NGS Pipeline

Next-generation sequencing (NGS) data analysis pipelines for RNA-seq, CUT&Tag, and CRISPR screening.

## Installation

### 1. Clone Repository
```bash
git clone https://github.com/soulong/NGS_pipeline.git
cd ngs_pipeline
```

### 2. Create Conda Environment
```bash
conda env create -f environment.yml
conda activate ngs
```

### 3. Verify Installation
```bash
# Check core tools
salmon --version
star --version
bowtie2 --version
```

## Quick Start

### RNA-seq Analysis
```bash
# Navigate to your dataset
cd /path/to/dataset

# Run pipeline
bash /path/to/ngs_pipeline/RNAseq/RNAseq_Step_1_Run.sh config.yml
```

### CUT&Tag Analysis
```bash
# Generate samplesheet
bash /path/to/ngs_pipeline/CutTag/run_cuttag.sh config.yml

# Run pipeline
bash /path/to/ngs_pipeline/CutTag/run_cuttag.sh config.yml
```

### CRISPR Screening
```bash
# Run MAGeCK analysis
bash /path/to/ngs_pipeline/CRISPR/MAGeCK_run.sh
```

## Pipelines

### RNA-seq Pipeline

**Workflow:**
```
FASTQ → QC → Quantification → Count Table → DE Analysis → Enrichment
```

**Steps:**
1. `Step 1:` QC & Quantification (Salmon/STAR)
2. `Step 2:` Gene/Transcript counting
3. `Step 3:` Differential expression (DESeq2)
4. `Step 4:` Functional enrichment analysis

**Input:**
- Raw FASTQ files
- `config.yml` - configuration file
- `samplesheet.csv` - sample metadata

**Output:**
- Count tables (TPM, raw counts)
- Differential expression results
- Enrichment analysis results

---

### CUT&Tag Pipeline

**Workflow:**
```
FASTQ → QC → Alignment → Peak Calling → Normalization → Differential Analysis
```

**Features:**
- Spike-free normalization
- Peak quantification
- Coverage visualization
- Differential binding analysis

**Input:**
- Raw FASTQ files
- `config.yml` - configuration file
- `samplesheet.csv` - sample metadata

**Output:**
- Aligned BAM files
- Peak files (BED/PEAK)
- Normalized count tables
- Quality control plots

---

### CRISPR Screening Pipeline

**Workflow:**
```
FASTQ → MAGeCK Count → MAGeCK Test → Gene Correction → Dependency Analysis
```

**Features:**
- MAGeCK analysis
- Gene symbol correction
- Gene dependency profiling
- Results visualization

**Input:**
- Raw FASTQ files from CRISPR screen
- `MAGeCK_mle_design.txt` - experimental design
- `gene_symbol_correction.xlsx` - correction table

**Output:**
- sgRNA count tables
- Gene-level statistics
- Gene dependency profiles

---

## File Structure

```
ngs_pipeline/
├── RNAseq/           # RNA-seq analysis
├── CutTag/           # CUT&Tag analysis
├── CRISPR/           # CRISPR screening
├── helper.sh         # Common utilities
├── make_index.sh     # Genome indexing
└── environment.yml   # Conda environment
```

## Requirements

- **OS:** Linux/macOS
- **RAM:** ≥16GB (≥64GB recommended for STAR)
- **Storage:** ≥50GB per RNA-seq sample

## License

MIT License - See [LICENSE](LICENSE) file for details.
