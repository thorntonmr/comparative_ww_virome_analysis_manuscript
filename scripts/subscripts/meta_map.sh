#!/bin/bash

set -euo pipefail

run_mapping="T" # only change for debugging / dev

# Create log file
log_filename="log-$(basename $0 .sh)_$(date +"%Y%m%d_%H%M%S").log"
exec > >(tee -a "$log_filename") 2>&1

echo "----- Starting script: $0 -----"
echo "Start time: $(date +"%Y"."%m"."%d %H":"%M":"%S")"
echo "   " 

# DEFINE VARIABLES
fastq_path="$1"
reference_genome="$2"
output_dir="$3"
strictness="$4"

if [[ "$strictness" == "--sensitive" ]]; then
    ANI_filter="90"
elif [[ "$strictness" == "--specific" ]]; then
    ANI_filter="95"
else
    echo "ERROR: Invalid strictness value."
    echo "  - Options are --sensitive | --specific"
    exit 1
fi

# Log settings
echo "----- Run settings -----"
echo "working_dir:              $(pwd)"
echo "reference_genome:         ${reference_genome}"
echo "fastq_path:               ${fastq_path}"
echo "output_dir:               ${output_dir}"
echo "Using strictness option:  $strictness"
echo "  - ANI_filter:  ${ANI_filter}"
echo "   "

mkdir -p ${output_dir}
mkdir -p ${output_dir}/mapping_stats

echo "----- Running Minimap2 Mapping -----"
if [[ "$run_mapping" = "T" ]]; then
        for fastq_file in "${fastq_path}"/*R1.fastq.gz; do
                samp=$(basename ${fastq_file} _R1.fastq.gz)
                samp=$(basename ${samp} ${fastq_path})

                echo "   "
                echo "Mapping ${samp}"
                echo "Time: $(date +"%Y"."%m"."%d %H":"%M":"%S")"

                # Check fastqs
                if [[ -s ${fastq_file} ]] && [[ -s ${fastq_path}/${samp}_R2.fastq.gz ]]; then
                        echo "  - Found ${samp} fastq"
                        echo "  - Mapping to ${reference_genome}"
                else
                        echo "  Tried to process: "
                        echo "     - R1: ${fastq_file}"
                        echo "     - R2: ${fastq_path}/${samp}_R2.fastq.gz"
                        echo "  - ERROR: ${samp} read file not found or empty. Exiting."
                        exit 1
                fi      

                # Check index (checks last file of index creation)
                if [[ -s ${reference_genome} ]]; then
                        # Map
                        minimap2 \
                                --secondary=yes \
                                -N 10 \
                                -p 0.6 \
                                --eqx \
                                -ax sr \
                                -o ${output_dir}/${samp}.sam \
                                ${reference_genome} \
                                ${fastq_file} \
                                ${fastq_path}/${samp}_R2.fastq.gz 
                else
                        echo "  ${reference_genome} not found. Exiting."
                        exit 1
                fi

                # Create sorted, indexed bams
                if [[ -s ${output_dir}/${samp}.sam ]]; then
                        echo "  - Converting sam to bam..."
                        samtools view -b -h ${output_dir}/${samp}.sam > "${output_dir}/${samp}_temp.bam" 
                        
                        echo "  - Setting fuzzy multimappers..."
                        perl ./subscripts/helper_scripts/set_fuzzy_multimappers.pl \
                                -l 1 \
                                -i "${output_dir}/${samp}_temp.bam" \
                                -st='--threads=4' |\
                        samtools view -b -h - > "${output_dir}/${samp}.bam"

                        echo "  - Filtering with msamtools filter..."
                        msamtools filter -h -b -l 60 -p "$ANI_filter" -z 80 "${output_dir}/${samp}.bam" | \
                        samtools sort -@4 -m4g > "${output_dir}/${samp}_mapped.bam"

                        echo "  - Indexing mapped bam..."
                        samtools index ${output_dir}/${samp}_mapped.bam
                else
                        echo "  WARNING: sam file ${samp}.sam does not exist or is empty. Continuing."
                fi

                # Clean up
                if [ -s ${output_dir}/${samp}.bam ]; then
                        echo "  - bam file successfully created and sorted. Removing sam and temp bam precusor."
                        rm "${output_dir}/${samp}.sam"
                        rm "${output_dir}/${samp}_temp.bam"
                else
                        echo "  WARNING: bam file creation error for ${samp}. Saving sam."
                fi

        done # close per fastq sample loop
else
        echo "WARNING: Mapping will be skipped. Continuning..."
fi
echo "  "

echo "----- Mapping complete.  Sorted bams are in ${output_dir}. -----"
