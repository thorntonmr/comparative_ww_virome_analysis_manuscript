#!/bin/bash

set -eu

# Check arguments
if [ "$1" != "bam" ] && [ "$1" != "fastq" ] && [ "$1" != "--help" ]; then
    echo "ERROR: Must specify bam | fastq | --help in script call. Exiting..."
    exit 1
fi

# Show usage
if [[ "$1" == "--help" ]]; then
    echo "This script runs QC either from bams (converts to fastq.gz), or fastq.gz files"
    echo "  - Two arguments are required. bam or fastq in the first position, and the directory path in the second position."
    echo "  - Fastq files should be in the *_R1.fastq.gz *_R2.fastq.gz format."
    exit 0
fi

# Create log file
log_filename="log-$(basename $0 .sh)_$(date +"%Y%m%d_%H%M%S").log"
exec > >(tee -a "$log_filename") 2>&1

echo "----- Starting script: $0 -----"
echo "Start time: $(date +"%Y"."%m"."%d %H":"%M":"%S")"
echo "   "

# Log settings
echo "----- Run Settings -----"
echo "Working directory = $(pwd)"
if [[ "$1" == "bam" ]]; then 
    echo "bam_dir = $2"
elif [[ "$1" == "fastq" ]]; then
    echo "fastq_dir = $2"
fi
echo "Run options: $1"
echo "  - Starting file type = $1"
echo "  - File directory = $2"
echo "   "

echo "Checking input files..."
# Check bam inputs
if [[ "$1" ==  "bam" ]]; then
    bam_dir="$2"
    if [[ ! -d "$bam_dir" ]]; then 
        echo "bam_dir ${bam_dir} does not exist.  Exiting..."
        exit 1
    elif [[ -z "$(find "$bam_dir" -mindepth 1 -maxdepth 1 -type f -name *.bam)" ]]; then
        echo "bam_dir ${bam_dir} is empty. Exiting..."
        exit 1
    else
        for i in "$bam_dir"/*.bam; do
            samp=$(basename "$i" "$bam_dir")
            echo "  - ${samp} found and ready for QC"
        done
    fi
fi
echo "   "

# Check fastq inputs
if [[ "$1" ==  "fastq" ]]; then
    fastq_dir="$2"
    if [[ ! -d "$fastq_dir" ]]; then 
        echo "fastq_dir ${fastq_dir} does not exist.  Exiting..."
        exit 1
    elif [[ -z "$(find "$fastq_dir" -mindepth 1 -maxdepth 1 -type f -name "*.fastq.gz")" ]]; then
        echo "Fastq dir is empty or file formats are wrong (should be R1/R2.fastq.gz). Exiting..."
        exit 1
    else
        for i in "$fastq_dir"/*R1.fastq.gz; do
            samp=$(basename $i _R1.fastq.gz)
            samp=$(basename $samp $fastq_dir)
            echo "  - $samp found and ready for QC..."
        done
    fi

fi

if [[ "$1" == "bam" ]]; then
    echo "Converting bams to paired end fastqs..."
    mkdir -p fastq
    for i in "$bam_dir"/*.bam; do
        samp=${i#*#}
        samp=$(basename $samp .bam)
        echo "  - Converting to fastq for: $samp"
        samtools fastq $i -1 fastq/${samp}_R1.fastq.gz -2 fastq/${samp}_R2.fastq.gz
        echo "   "
    done

    # Set fastq dir
    fastq_dir="fastq"
else
    echo "Starting from fastq.gz files. No need to convert..."
    echo "   "
fi

mkdir -p trimmed_fastq
mkdir -p fastp_reports

echo "Running Fastp QC..."
for i in ${fastq_dir}/*R1.fastq.gz; do
    samp=$(basename $i _R1.fastq.gz)

    echo "  - Running fastp for ${samp}..."
    fastp \
        --low_complexity_filter \
        --correction \
        --complexity_threshold 30 \
        --length_required 50 \
        --average_qual 25 \
        --qualified_quality_phred 25 \
        --unqualified_percent_limit 60 \
        --cut_mean_quality 25 \
        --cut_front \
        --cut_tail \
        --dedup \
        --dup_calc_accuracy 5 \
        --in1 ${i} \
        --in2 ${fastq_dir}/${samp}_R2.fastq.gz \
        --out1 trimmed_fastq/${samp}_trimmed_R1.fastq.gz \
        --out2 trimmed_fastq/${samp}_trimmed_R2.fastq.gz \
        --html fastp_reports/${samp}_fastp.html \
        --json fastp_reports/${samp}_fastp.json
    echo "  "
    
    # Check if reads remain after filtering
    if [[ -s trimmed_fastq/${samp}_trimmed_R1.fastq.gz ]] && [[ -s trimmed_fastq/${samp}_trimmed_R2.fastq.gz ]]; then
        # Run qc on pre and post filtered reads
        echo "  - Running fastqc..."
        fastqc \
            ${i} \
            ${fastq_dir}/${samp}_R2.fastq.gz \
            trimmed_fastq/${samp}_trimmed_R1.fastq.gz \
            trimmed_fastq/${samp}_trimmed_R2.fastq.gz \
            --outdir fastp_reports \
            --threads 8 \
            --quiet
    else
        echo "  - WARNING: no reads remain after filtering for $samp. Removing empty file..."
        rm trimmed_fastq/${samp}_trimmed_R1.fastq.gz
        rm trimmed_fastq/${samp}_trimmed_R2.fastq.gz
    fi
done

echo "Consolidating QC reports with MultiQC"
multiqc fastp_reports -o fastp_reports

echo "----- $0 finished successfully -----"