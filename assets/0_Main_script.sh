#!/bin/bash
CONDA_BASE_PATH=$(conda info --base)
source $CONDA_BASE_PATH/etc/profile.d/conda.sh

#####################################
# The main script will run on the frontend (it is not consuming many resources).
# All processes will run on different nodes.
####################################


# paths
NF_FOLDER=/mnt/ngs_ricerca/NEXTFLOW/
UsefulData=/mnt/ngs_ricerca/Software/
work_dir=/mnt/ngs_ricerca/NEXTFLOW/nextflow_temp/
conf_file=/mnt/ngs_ricerca/NEXTFLOW/local.config

######## PARAMS ########
sample_file=/mnt/ngs_ricerca/dichiarop/Master_scripts/LMD/Sample.txt
outdir=/mnt/ngs_ricerca/dichiarop/Data/LMD_data/OUT/
seq_mode=PE
genome_fasta=$UsefulData/reference_genome/GRCh38.primary_assembly.genome.fa
transcriptome_fasta=$UsefulData/reference_transcriptome/Gencode_v47_GRCh38_p14/gencode.v47.transcripts.fa
index=$UsefulData/reference_genome/kallisto_index/kallisto.v0.48_gencode.v47_hg38_k31/kallisto
GTF=$UsefulData/reference_genome/gencode.v47.basic.annotation.gtf   #gencode.v47.primary_assembly.annotation.gtf
reference=$UsefulData/reference_genome/hg38_annotation/HG38_merged_annotation.txt
########################

# Where will be stored work files, it can be deleted once the script is finished
mkdir -p $work_dir

project_name=test
work_dir_project=$work_dir/$project_name/

mkdir -p $work_dir_project
cd $work_dir_project


# -------------------------
# Run Nextflow workflow
# -------------------------
### NOTE ###
#--email pierluigi.dichiaro@ausl.re.it \   #if MTA (Mail Transfer Agent) is installed
#--genome GRCh38  \  #GRCh38 --> NCBI  #hg38 --> UCSC   ## if do not you have genome .fasta (not raccomended)
conda activate nextflow

export NXF_ASSETS=/mnt/ngs_ricerca/NEXTFLOW/Nextflow_pipeline

NXF_VER=25.04.7 nextflow run pdichiaro/LMDseq -r main \
            --input $sample_file \
            --outdir $outdir \
            --fasta $genome_fasta \
            --transcript_fasta $transcriptome_fasta \
            --index $index \
            --seq_mode $seq_mode \
            --gtf $GTF \
            --reference $reference \
            --gencode True \
            --trimmer trimgalore \
            --pseudo_aligner kallisto \
            --skip_gtf_filter True \
            --skip_fastqc False \
            --skip_trimming False \
            --skip_pseudo_alignment False \
            --save_reference False \
            -profile singularity \
            -c $conf_file \
            -process.echo \
            -resume \
            -with-report report.html \
            -with-trace trace.txt \
            -with-timeline timeline.html \
            -with-dag flowchart.html


if [ $? -eq 0 ]; then
    echo "Nextflow finished successfully"
    scp -r {./report*,./timeline*,./trace*,./flowchart*} $outdir
    rm -rf $work_dir_project
else
    echo "Nextflow encountered an error"
fi

conda deactivate


