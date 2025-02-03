#!/bin/bash

set -eu

# Create log file
log_filename="log-$(basename $0 .sh)_$(date +"%Y%m%d_%H%M%S").log"
exec > >(tee -a "$log_filename") 2>&1

echo "----- Starting script: $0 -----"
echo "Start time: $(date +"%Y"."%m"."%d %H":"%M":"%S")"
echo "   " 

# Set variables
fastq_dir="host-filt_fastq"
subset_fasta_dir="subset_fastas"
helper_scripts="subscripts/helper_scripts"
strictness="$1"

map="T" # only change for rerunning stat section and/or debugging

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
echo "fastq_dir = ${fastq_dir}"
echo "subset_fasta_dir = ${subset_fasta_dir}"
echo "Using strictness option = ${strictness}"
echo "  - ANI_filter = ${ANI_filter}"
echo "   "

# Check fastq inputs
if [[ -z "$(find "$fastq_dir" -mindepth 1 -maxdepth 1 -type f -name "*fastq.gz")" ]]; then
    echo "ERROR: ${fastq_dir} not found or empty. Exiting..."
    echo "  - Looking for: $(find "$fastq_dir" -mindepth 1 -maxdepth 1 -type f -name "*fastq.gz")"
    exit 1
else
    echo "Found fastq directory..."
    echo "  - $(find ${subset_fasta_dir} -name "*subset.fasta" | wc -l) samples passing filter ready for processing."
fi

# Remap to subset reference
mkdir -p filtered_mapping_results

echo "Mapping to filtered references..."
if [[ "$map" = "T" ]]; then
    for fasta in "${subset_fasta_dir}"/*_subset.fasta; do
        samp=$(basename ${fasta} _mapped_fasta_subset.fasta)
        samp=$(basename ${samp} ${subset_fasta_dir})

        # Check fastqs
        if [[ -s "${fastq_dir}"/${samp}_R1.fastq.gz ]] && [[ -s "${fastq_dir}"/${samp}_R2.fastq.gz ]]; then
            echo "Found ${samp} fastq..." 
            echo "Mapping to $fasta"
        else
            echo "  ERROR: ${samp} read file not found or empty. Exiting."            
            echo "  Tried to process: "
            echo "     - R1: ${fastq_dir}/${samp}_R1.fastq.gz"
            echo "     - R2: ${fastq_dir}/${samp}_R2.fastq.gz"
            exit 1
        fi  

        # Check index (checks last file of index creation)
        if [[ -s "${subset_fasta_dir}/${samp}_mapped_fasta_subset.fasta" ]]; then
            # Map
            minimap2 \
                    --secondary=yes \
                    --eqx \
                    -ax sr \
                    "$fasta" \
                    "${fastq_dir}"/${samp}_R1.fastq.gz \
                    "${fastq_dir}"/${samp}_R2.fastq.gz \
                    -o filtered_mapping_results/${samp}_tmp.sam

            echo "  - Sorting reads by name..."
            samtools view -b -h filtered_mapping_results/${samp}_tmp.sam | samtools sort -n -@6 -m4g > filtered_mapping_results/${samp}_tmp.bam

            echo "  - Filtering reads + sorting by position..."
            msamtools filter -h -b -l 60 -p "$ANI_filter" -z 80 --besthit filtered_mapping_results/${samp}_tmp.bam | \
            samtools sort -@6 -m4g > filtered_mapping_results/${samp}_mapped.bam

            echo "  - Indexing sorted bam..."
            samtools index filtered_mapping_results/${samp}_mapped.bam
            
            # Check success
            if [[ -s "filtered_mapping_results/${samp}_mapped.bam" ]]; then
                echo "  - Removing temp files..."
                rm filtered_mapping_results/*tmp*
                echo "   "
            else
                echo "  WARNING: filtered_mapping_results/${samp}_mapped.bam not found. Continuing..."
                echo "   "
            fi
        else
            echo "  $fasta not found. Exiting."
            exit 1
        fi
    done
else
    echo "WARNING: Skipping mapping"
fi

# Get stats
echo "----- Mapping finished, getting stats -----"
mkdir -p filtered_mapping_results/mapping_stats/individual_files
mkdir -p filtered_mapping_results/mapping_stats

# Sequencing stats
for i in filtered_mapping_results/*_mapped.bam; do
    samp=$(basename $i .bam)

    echo "Processing stats for ${samp}"

    # Reads per library
    echo "  Counting reads per library..."
    echo -e "${samp}\t$(samtools view -c "${i}")" >> "filtered_mapping_results/mapping_stats/reads_per_library.txt"
    
    # Per base coverage
    echo "  Calculating read depth..."
    samtools depth "$i" > "filtered_mapping_results/mapping_stats/${samp}_coverage.tsv"

    # Compile number of non-primary mappers per accession
    echo "  Counting secondary mappers and alignment stats..."
    Rscript "${helper_scripts}/get_secondary_mappers.R" filtered_mapping_results ${samp}
    echo "   "

done # close per sample loop

echo "   "
echo "----- Finished gathering stats, running summarize_stats.R -----"

# Run summary script. Should be kept in standard directory.
Rscript "$helper_scripts"/summarize_stats.R filtered_mapping_results/mapping_stats

# Move individual files to sub-directory and keep summary files separate
if [[ -s filtered_mapping_results/mapping_stats/stat_summary.csv ]]; then
        mkdir -p final_stats
        mv filtered_mapping_results/mapping_stats/*coverage.tsv filtered_mapping_results/mapping_stats/individual_files
        mv filtered_mapping_results/mapping_stats/*idxstat.tsv filtered_mapping_results/mapping_stats/individual_files
        mv filtered_mapping_results/mapping_stats/stat_summary.csv ./final_stats/final_mapping_stats.csv
else
        echo "ERROR: summarizing statistics failed. Exiting."
        exit 1
fi

echo "   "
echo "----- Stats compiled,  $0 finished succesfully. -----"
echo 'End time:' `date +'%Y-%m-%d %T'`