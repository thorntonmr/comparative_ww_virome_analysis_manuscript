#!/bin/bash

set -eu

# Create log file
log_filename="log-$(basename $0 .sh)_$(date +"%Y%m%d_%H%M%S").log"
exec > >(tee -a "$log_filename") 2>&1

echo "----- Starting script: $0 -----"
echo "Start time: $(date +"%Y"."%m"."%d %H":"%M":"%S")"
echo "   " 

# Define variables
output_dir="$1"
helper_scripts="subscripts/helper_scripts"

mkdir -p ${output_dir}/mapping_stats/individual_files
mkdir -p ${output_dir}/bams_unfiltered
mkdir -p ${output_dir}/bams_mapped

# Sequencing stats
for i in ${output_dir}/*mapped.bam; do
    samp=$(basename $i .bam)

    echo "Processing stats for ${samp}"

    # Reads per library
    echo "  - Counting reads per library..."
    echo -e "${samp}\t$(samtools view -c "$(basename $i _mapped.bam).bam")" >> "${output_dir}/mapping_stats/reads_per_library.txt"
    
    # Per base coverage
    echo "  - Calculating read depth..."
    samtools depth ${output_dir}/${samp}.bam > ${output_dir}/mapping_stats/${samp}_coverage.tsv

    # Compile number of non-primary mappers per accession
    echo "  - Counting secondary mappers and alignment stats..."
    Rscript "$helper_scripts"/get_secondary_mappers.R ${output_dir} ${samp}
    echo "   "

done # close per sample loop

# Clean up empty files
for file in "${output_dir}"/mapping_stats/*.tsv; do
    if [[ ! -s $file ]]; then
        echo "  - WARNING: removing empty file $file" 
        rm $file
    fi
done

echo "----- Finished gathering stats, running summarize_stats.R -----"

# Run summary script. Should be kept in standard directory.
Rscript "$helper_scripts"/summarize_stats.R ${output_dir}/mapping_stats

# Move individual files to sub-directory and keep summary files separate
if [[ -s ${output_dir}/mapping_stats/stat_summary.csv ]]; then
        mv ${output_dir}/mapping_stats/*coverage.tsv ${output_dir}/mapping_stats/individual_files
        mv ${output_dir}/mapping_stats/*idxstat.tsv ${output_dir}/mapping_stats/individual_files
else
        echo "ERROR: summarizing statistics failed. Exiting."
        exit 1
fi

# Consolidate bams in sub-directories
mv ${output_dir}/*_mapped.bam* ${output_dir}/bams_mapped
mv ${output_dir}/*.bam* ${output_dir}/bams_unfiltered

echo "   "
echo "----- Stats compiled, script finished succesfully. -----"
echo 'End time:' `date +'%Y-%m-%d %T'`
