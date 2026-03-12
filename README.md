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

## Requirements

- **OS:** Linux/macOS
- **RAM:** ≥16GB (≥64GB recommended for STAR)
- **Storage:** ≥50GB per RNA-seq sample

## License

MIT License - See [LICENSE](LICENSE) file for details.
