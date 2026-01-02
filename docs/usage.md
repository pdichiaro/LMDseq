# LMDseq Usage Guide

## Overview

The **LMDseq** pipeline is a specialized bioinformatics workflow designed for RNA sequencing analysis of **low input samples obtained through Laser Microdissection (LMD)**. This pipeline addresses the unique analytical challenges of processing RNA-seq data from microdissected tissue samples, which typically yield very small amounts of RNA and require optimized computational approaches.

## Pipeline Features

- **LMD-Optimized**: Specifically designed for low input RNA samples from laser microdissection
- **Kallisto Pseudo-alignment**: Fast and accurate quantification optimized for LMD samples  
- **Quality Control**: Comprehensive FastQC analysis of raw and trimmed reads
- **Adapter Trimming**: TrimGalore with LMD-specific quality parameters
- **Flexible Output**: Minimal (gene_id + counts) or enriched (full annotations) matrices
- **Custom References**: Support for custom genome, transcript, and annotation files
- **DESeq2 QC**: Size factor normalization, VST, PCA, and correlation analysis
- **Reproducible**: Containerized environment ensuring consistent results
- **Streamlined**: No UMI processing or rRNA depletion (not needed for LMD samples)


## Input Preparation

### Sample.txt File Format

The pipeline requires a **tab-separated** sample sheet (`Sample.txt`) located in the `assets/` directory. This file contains metadata about your LMD samples and the location of FASTQ files.

#### Required Columns:

| Column | Description | Example | Notes |
|--------|-------------|---------|-------|
| `SID` | Sample ID (unique identifier) | LCM_1, LCM_2 | Must be unique across all samples |
| `SHORT_NAME` | Short sample name for output files | 20942_1D3_A1 | Used in final sample naming |
| `REPLICATE` | Replicate identifier | R1, R2, R3 | Biological or technical replicate |
| `RID` | Run ID from sequencing | S61303 | Sequencing run identifier |
| `strandedness` | Library strandedness | reverse, forward, unstranded | Must match library prep protocol |
| `FASTQ_FOLDER` | Full path to directory containing FASTQ files | /path/to/fastq/files | Must be absolute path |

#### Example Sample.txt:
```
SID	SHORT_NAME	REPLICATE	RID	strandedness	FASTQ_FOLDER
LCM_1	20942_1D3_A1	R1	S61303	reverse	/data/fastq/Sample_S61303_LCM_1
LCM_2	20942_1D3_B1	R2	S61304	reverse	/data/fastq/Sample_S61304_LCM_2
LCM_3	20942_1D3_C1	R3	S61305	reverse	/data/fastq/Sample_S61305_LCM_3
```

### Creating Your Sample.txt File

1. **Prepare your data structure:**
   - Organize FASTQ files in separate directories for each sample
   - Ensure FASTQ files follow naming convention: `*_R1_*.fastq.gz` and `*_R2_*.fastq.gz` for paired-end data
   - Use consistent naming across all samples

2. **Create the Sample.txt file:**
```bash
# Navigate to the assets directory
cd assets/

# Create or edit Sample.txt
nano Sample.txt
```

3. **Fill in the required information:**
   - Use **absolute paths** for `FASTQ_FOLDER`
   - Ensure `strandedness` matches your library preparation protocol
   - Use consistent naming conventions for `SHORT_NAME`
   - Verify all FASTQ directories exist and are accessible

### FASTQ File Organization

Your FASTQ files should be organized as follows:
```
/path/to/data/
├── Sample_S61303_LCM_1/
│   ├── sample1_R1_001.fastq.gz
│   └── sample1_R2_001.fastq.gz
├── Sample_S61304_LCM_2/
│   ├── sample2_R1_001.fastq.gz
│   └── sample2_R2_001.fastq.gz
└── Sample_S61305_LCM_3/
    ├── sample3_R1_001.fastq.gz
    └── sample3_R2_001.fastq.gz
```

## Running the Pipeline

The LMDseq pipeline offers flexible reference genome options and output formats to suit different analysis needs.

### Quick Start Examples

#### **1. Basic Analysis (Minimal Output)**
```bash
nextflow run pdichiaro/LMDseq \
  --input assets/Sample.txt \
  --genome GRCh38 \
  --outdir results \
  -profile docker
```
- **Output**: Clean gene expression matrix (gene_id + counts)
- **Use case**: Quick quantification and exploratory analysis

#### **2. Custom Genome Files (Minimal Output)**
```bash
nextflow run pdichiaro/LMDseq \
  --input assets/Sample.txt \
  --fasta /path/to/genome.fa \
  --gtf /path/to/annotation.gtf \
  --outdir results \
  -profile docker
```
- **Output**: Minimal matrix with custom genome
- **Use case**: Non-standard organisms or custom assemblies

#### **3. Enriched Analysis (Full Annotations)**
```bash
nextflow run pdichiaro/LMDseq \
  --input assets/Sample.txt \
  --fasta /path/to/genome.fa \
  --gtf /path/to/annotation.gtf \
  --reference /path/to/gene_annotations.txt \
  --outdir results \
  -profile docker
```
- **Output**: Full gene matrix with symbols, names, coordinates
- **Use case**: Publication-ready analysis with detailed annotations

#### **4. Transcript-Only Mode**
```bash
nextflow run pdichiaro/LMDseq \
  --input assets/Sample.txt \
  --transcript_fasta /path/to/transcripts.fa \
  --outdir results \
  -profile docker
```
- **Output**: Direct transcript quantification
- **Use case**: Custom transcriptomes or pre-built indices

### Reference Genome Options

| Parameter | Type | Description | Example |
|-----------|------|-------------|---------|
| `--genome` | String | iGenomes reference ID | `GRCh38`, `GRCm39`, `R64-1-1` |
| `--fasta` | Path | Custom genome FASTA file | `/path/to/genome.fa` |
| `--gtf` | Path | Gene annotation GTF file | `/path/to/genes.gtf` |
| `--transcript_fasta` | Path | Custom transcript FASTA file | `/path/to/transcripts.fa` |
| `--reference` | Path | Gene annotation reference file **(recommended)** | `/path/to/annotations.txt` |
| `--index` | Path | Pre-built Kallisto index | `/path/to/kallisto.idx` |

### Core Pipeline Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--input` | Path | **Required** | Sample metadata file (Sample.txt) |
| `--outdir` | Path | **Required** | Output directory path |
| `--pseudo_aligner` | String | `kallisto` | Pseudo-alignment tool (only kallisto supported) |
| `--skip_pseudo_alignment` | Boolean | `false` | Skip pseudo-alignment step |

###  Kallisto-Specific Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--extra_kallisto_quant_args` | String | `'-b 50 --single-overhang --genomebam'` | Additional Kallisto quantification arguments |
| `--extra_kallisto_index_args` | String | `null` | Additional Kallisto indexing arguments |
| `--kallisto_quant_fraglen` | Integer | `200` | Fragment length for single-end mode |
| `--kallisto_quant_fraglen_sd` | Integer | `200` | Standard deviation for fragment length |
| `--kallisto_kmer_size` | Integer | `31` | K-mer size for Kallisto index |

###  Quality Control & Trimming Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--trimmer` | String | `trimgalore` | Trimming tool (only trimgalore supported) |
| `--extra_trimgalore_args` | String | `'--clip_R2 3 --quality 20 --stringency 3 --length 20'` | TrimGalore arguments |
| `--skip_fastqc` | Boolean | `false` | Skip FastQC analysis |
| `--skip_trimming` | Boolean | `false` | Skip adapter trimming |

###  DESeq2 QC Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--min_reads` | Integer | `10` | Minimum read count threshold for gene filtering |

###  Reference Genome Processing

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--gencode` | Boolean | `false` | Use GENCODE transcript processing |
| `--skip_gtf_filter` | Boolean | `true` | Skip GTF filtering step |
| `--save_reference` | Boolean | `false` | Save reference files to output directory |
| `--igenomes_ignore` | Boolean | `false` | Ignore iGenomes reference configuration |

###  Advanced Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--seq_mode` | String | `SE` | Sequencing mode (SE=single-end, PE=paired-end) |
| `--bin_size` | Integer | `1` | Bin size for coverage analysis |
| `--chromosomes` | Path | `null` | Chromosome list file |


### Output Formats

The pipeline generates different output formats depending on whether a reference annotation file is provided:

### **Minimal Output (No --reference)**
```
gene_id        SAMPLE1  SAMPLE2  SAMPLE3
ENSG00000001   245      189      312
ENSG00000002   156      203      89
ENSG00000003   89       134      201
```
### **Enriched Output (With --reference)**
```
gene_id        Symbol  Name                  type_of_gene    seqnames  start     SAMPLE1  SAMPLE2  SAMPLE3
ENSG00000001   TP53    tumor protein p53     protein_coding  chr17     7565097   245      189      312  
ENSG00000002   BRCA1   breast cancer gene 1  protein_coding  chr17     43044295  156      203      89
ENSG00000003   TERT    telomerase reverse    protein_coding  chr5      1253147   89       134      201
```

| Parameter | Description | Required | Default |
|-----------|-------------|----------|---------|
| `--genome` | iGenomes genome ID (e.g., GRCh38, GRCm39) | ❌ | `null` |
| `--fasta` | Custom genome FASTA file | ❌ | `null` |
| `--gtf` | Gene annotation GTF file | ❌ | `null` |
| `--transcript_fasta` | Custom transcript FASTA file | ❌ | `null` |
| `--reference` | Gene annotation reference file (optional) | ❌ | `null` |
| `--index` | Pre-built Kallisto index | ❌ | `null` |

### Core Parameters

| Parameter | Description | Required | Default |
|-----------|-------------|----------|---------|
| `--input` | Sample metadata file (Sample.txt) | ✅ | - |
| `--outdir` | Output directory | ✅ | - |
| `--pseudo_aligner` | Pseudo-alignment tool | ❌ | `kallisto` |
| `--seq_mode` | Sequencing mode (SE/PE) | ❌ | `SE` |
| `--min_reads` | Minimum read threshold for gene filtering | ❌ | `10` |

### Kallisto-Specific Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--kallisto_quant_fraglen` | Fragment length for single-end reads | `200` |
| `--kallisto_quant_fraglen_sd` | Fragment length standard deviation | `200` |
| `--kallisto_kmer_size` | K-mer size for indexing | `31` |
| `--extra_kallisto_quant_args` | Additional Kallisto quantification arguments | `'-b 50 --single-overhang --genomebam'` |
| `--extra_kallisto_index_args` | Additional Kallisto indexing arguments | `null` |
```

### Test Run
```bash
nextflow run pdichiaro/LMDseq \
  -profile test_lmd \
  --outdir test_results
```

## Parameters

### Essential Parameters

#### Input/Output
- `--input`: Path to Sample.txt file (required)
- `--outdir`: Output directory for results (required)

#### Reference Genome
- `--genome`: iGenomes reference (e.g., GRCh38, GRCm38, GRCz10)
- `--fasta`: Custom genome FASTA file
- `--gtf`: Custom GTF annotation file
- `--kallisto_index`: Pre-built Kallisto index
- `--reference`: Reference annotation file for enhanced gene expression matrix

### LMD-Specific Parameters

#### Kallisto Quantification
- `--aligner`: Alignment method (default: `kallisto`)
- `--extra_kallisto_quant_args`: Additional Kallisto arguments (default: `-b 50 --single-overhang --genomebam`)
- `--extra_kallisto_index_args`: Additional Kallisto index arguments
- `--kallisto_quant_fraglen`: Fragment length for single-end reads (default: 200)
- `--kallisto_quant_fraglen_sd`: Fragment length standard deviation (default: 200)
- `--kallisto_kmer_size`: K-mer size for Kallisto index (default: 31)

#### DESeq2 Analysis
- `--min_reads`: Minimum read count threshold for gene filtering (default: 10)

#### Read Processing
- `--trimmer`: Trimming tool (default: `trimgalore`)
- `--extra_trimgalore_args`: Additional TrimGalore arguments (default: `--clip_R2 3 --quality 20 --stringency 3 --length 20`)

### Process Control

#### Quality Control
- `--skip_fastqc`: Skip FastQC analysis

#### Process Skipping
- `--skip_trimming`: Skip adapter trimming step
- `--skip_gtf_filter`: Skip GTF filtering (default: true)


## Output Structure

The pipeline generates a structured output directory with numbered folders for easy navigation:

```
results/
├── 1_FASTQ/                    # Merged FASTQ files
├── 2_QC/                       # FastQC quality control reports
├── 3_TRIM/                     # Adapter trimming outputs and reports
├── 4_BAM/                      # Kallisto quantification results and pseudo-BAM files
├── 5_Bw/                       # BigWig coverage files for genome browser visualization
├── 6_Norm_folder/              # Normalization plots and read distribution analysis
├── 7_Counts_folder/            # Gene expression matrices (raw and normalized)
├── 8_Quality_folder/           # Advanced QC analysis (PCA, correlations, distances)
└── pipeline_info/              # Pipeline execution information and resource usage
```

### Output Matrix Formats

The pipeline generates different matrix formats depending on whether a reference annotation file is provided:

### Key Output Files

#### **Gene Expression Matrices**
| File | Description | Format |
|------|-------------|--------|
| `EX_reads_RAW.txt` | Raw count matrix | Minimal or Enriched |
| `EX_reads_NORM_filt.txt` | DESeq2 size-factor normalized counts | Same as input |
| `rLog_reads.txt` | Variance stabilized (VST) data | Same as input |
| `Normalisation_Parameters.txt` | DESeq2 size factors | Sample-wise normalization factors |

#### **Quality Control Reports**
| File | Description |
|------|-------------|
| `PCA_rlogTransformed.pdf` | Principal component analysis plots |
| `Heatmap_sampleTosample_correlation.pdf` | Sample correlation matrix |
| `Heatmap_sampleTosample_distances.pdf` | Sample distance clustering |

#### **Coverage and Visualization**
| File | Description |
|------|-------------|
| `*.bw` | BigWig coverage files for genome browsers |
| `*.bam` | Kallisto pseudo-alignments (genomic coordinates) |

## Advanced Usage

### **Profile Options**
```bash
# Docker (recommended)
-profile docker

# Singularity
-profile singularity  

# Conda
-profile conda

# Test profiles
-profile test_lmd    # Small test dataset
```

### **Custom Configuration**
```bash
# Custom resource limits
nextflow run pdichiaro/LMDseq \
  --input samples.txt \
  --genome GRCh38 \
  --max_memory '32.GB' \
  --max_cpus 16 \
  --outdir results
```

### **Troubleshooting**

#### **Common Issues:**
1. **Missing FASTQ files**: Ensure `FASTQ_FOLDER` paths in Sample.txt are correct and accessible
2. **Memory errors**: Increase `--max_memory` for large genomes 
3. **No genes found**: Check GTF file format and genome compatibility
4. **Reference errors**: Ensure reference file matches your genome annotation

#### **Getting Help:**
```bash
# View all parameters
nextflow run pdichiaro/LMDseq --help

# Test pipeline
nextflow run pdichiaro/LMDseq -profile test_lmd --outdir test_results

# Check version
nextflow run pdichiaro/LMDseq --version
```


## Support and Citation

### Getting Help
- **Issues**: [GitHub Issues](https://github.com/pdichiaro/LMD_rnaseq/issues)
- **Documentation**: Complete documentation in `docs/` directory
- **Contact**: Pierluigi Di Chiaro (AUSL-IRCCS Reggio Emilia)


### Citation
If you use LMD_RNAseq in your research, please cite:

> **Di Chiaro P, et al.** "Mapping functional to morphological variation reveals the basis of regional extracellular matrix subversion and nerve invasion in pancreatic cancer." *Cancer Cell* 2024.
> **Di Chiaro P, et al.** "A framework to mine laser microdissection-based omics data and uncover regulators of pancreatic cancer heterogeneity." *Gigascience* 2025.
