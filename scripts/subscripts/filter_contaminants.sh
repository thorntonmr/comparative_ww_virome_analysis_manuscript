#!/bin/bash

set -eu
set -o pipefail

# Check if the correct number of arguments is provided
if [[ "$#" -ne 2 ]]; then
    echo "   "
    echo "ERROR: Invalid number of arguments"
    echo "   "
    exit 1
fi

# Check fastq_dir path
if [[ ! -d $1 ]]; then
    echo "   "
    echo "ERROR: fastq_dir does not exist."
    exit 1
fi

# Create log file
log_filename="log-$(basename $0 .sh)_$(date +"%Y%m%d_%H%M%S").log"
exec > >(tee -a "$log_filename") 2>&1

echo "Starting script: $0"
echo "Start time: $(date +"%Y"."%m"."%d %H":"%M":"%S")"

echo "   "
echo "----- Using run options -----"
echo "Working directory: $(pwd)"
echo "Fastq diretory: " $1
echo "Reference genome: " $2

# Set variables
fastq_dir=$1
bt2_idx=$2

# Check inputs
echo "----- Checking input files -----"
for fastq in ${fastq_dir}/*.fastq.gz; do
    if [[ -s $fastq ]]; then
        if [[ $fastq =~ R1 || $fastq =~ R2 ]]; then
            echo "  - $(basename $fastq .fastq.gz) ready for host removal"
        else
            echo "${fastq} exists but improper format. Exiting"
            exit 1
        fi
    else
        echo "${fastq} does not exist. Exiting"
        exit 1
    fi
done
echo "  "

# Create output dir
mkdir -p ./host-filt_fastq

echo "----- Performing paired end host removal -----"
# Remove host reads
for fastq in ${fastq_dir}/*R1.fastq.gz; do
    samp=$(basename $fastq _R1.fastq.gz)
    echo "   "
    echo "Filtering reads from ${samp}..."

    bowtie2 \
        -x ${bt2_idx} \
        --very-fast \
        -1 ${fastq} \
        -2 ${fastq_dir}/${samp}_R2.fastq.gz \
        --un-conc-gz ${samp}_host_filtered.fastq.gz > ${samp}_host_reads.sam

    # Check output
    if [[ -s ${samp}_host_filtered.fastq.2.gz ]]; then
        echo "Moving ${samp} to output directory."
        mv ${samp}_host_filtered.fastq* ./host-filt_fastq
        rm ${samp}_host_reads.sam
    else
        echo "WARNING: Filtering failed for ${samp}."
    fi
done

# Reformat output fastq filenames
for i in ./host-filt_fastq/*1.gz; do
    samp=$(basename $i .fastq.1.gz)
    new_file="${samp}_R1.fastq.gz"

    echo "  Renaming ${i} to ${new_file}"
    mv "$i" ./host-filt_fastq/"$new_file"
done

for i in ./host-filt_fastq/*2.gz; do
    samp=$(basename $i .fastq.2.gz)
    new_file="${samp}_R2.fastq.gz"

    echo "  Renaming ${i} to ${new_file}"
    mv "$i" ./host-filt_fastq/"$new_file"
done

echo "   "
echo "---- Host removal process completed -----"

if [[ ! -d ./fastp_reports ]]; then
    mkdir -p ./fastp_reports
fi

cp "$log_filename" ./fastp_reports
multiqc fastp_reports -o fastp_reports --force

echo 'End time:' `date +'%Y-%m-%d %T'`