#!/bin/bash

set -eu

if [ "$1" == "--help" ]; then
    echo "Salmon quantification script requires 2 arguments."
    echo "  - First argument: path to folder with mapped bam files"
    echo "  - Second argument: path to folder with subsetted fasta files"
    exit 0
elif [ "$#" -ne 2 ]; then
    echo "ERROR: Invalid number of arguments. Script requires 2 arguments."
    echo "  - First argument: path to folder with mapped bam files"
    echo "  - Second argument: path to folder with subsetted fasta files"
    exit 1
fi

# Create log file
log_filename="log-$(basename $0 .sh)_$(date +"%Y%m%d_%H%M%S").log"
exec > >(tee -a "$log_filename") 2>&1

echo "----- Starting script: $0 -----"
echo "Start time: $(date +"%Y"."%m"."%d %H":"%M":"%S")"
echo "   " 

# Set variables
bam_dir="$1"
fasta_subset_dir="$2"

# Log settings
echo "----- Run Settings -----"
echo "bam_dir = ${bam_dir}"
echo "fasta_subset_dir = ${fasta_subset_dir}"
echo "   "

# Check inputs
echo "Checking inputs..."
if [[ ! -d "$bam_dir" ]]; then 
    echo "  - ERROR: bam_dir ${bam_dir} does not exist.  Exiting..."
    exit 1
elif [[ -z "$(find "$bam_dir" -mindepth 1 -maxdepth 1 -type f -name *.bam)" ]]; then
    echo "  - ERROR: bam_dir ${bam_dir} is empty. Exiting..."
    exit 1
else
    for i in "$bam_dir"/*.bam; do
        samp=$(basename "$i" "$bam_dir")
        echo "  - $samp found and ready for quantification..."
    done
fi

# Covert bams to fastq
echo "Converting mapped bams back to fastq..."
mkdir -p ./temp_fastq
for i in "$bam_dir"/*.bam; do 
    samp="$(basename $i .bam)"
    echo "  Converting $samp to fastq..."
    samtools fastq $i -1 temp_fastq/${samp}_R1.fastq.gz -2 temp_fastq/${samp}_R2.fastq.gz
done

mkdir -p salmon_quants
mkdir -p salmon_quants/all_sample_quants

# Run Salmon Quantifications
echo "----- Starting salmon -----"
mkdir -p salmon_index_temp

for i in "$bam_dir"/*.bam; do
    samp=$(basename $i .bam)
    reference="${fasta_subset_dir}/${samp}_fasta_subset.fasta"
    
    echo "  Quantifying ${samp}..."
    echo "  Using reference ${reference}..."

    echo "  Preparing index..."
    salmon index -t "$reference" -i salmon_index_temp/${samp}_index -k 31 

    salmon quant \
        --index salmon_index_temp/"${samp}_index" \
        --libType A \
        --seqBias \
        --gcBias \
        -1 temp_fastq/${samp}_R1.fastq.gz \
        -2 temp_fastq/${samp}_R2.fastq.gz \
        --recoverOrphans \
        --softclipOverhangs \
        --meta \
        --writeMappings salmon_quants/${samp}_quants/${samp}_salmon_mapped.sam \
        --writeQualities \
        -o salmon_quants/${samp}_quants

    # Copy quant files to common folder
    cp salmon_quants/${samp}_quants/quant.sf salmon_quants/all_sample_quants/${samp}_quant.sf   
done

echo "----- Salmon quantification successfull -----"
echo "End time: $(date +"%Y"."%m"."%d %H":"%M":"%S")"