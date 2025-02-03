#!/bin/bash

set -euo pipefail

version="metamap_v1.3"

# Set Variables
starting_format="bam"                   # input file format (must be either .bam or _R*.fastq.gz)
starting_dir="path/to/starting/dir"     # directory path for raw input files (bam or fastq)
contaminant_dir="path/to/bt2_ref"       # path to bowtie2 index for contaminant reference
mapping_reference="path/to/ref/fasta"   # path to reference sequence(s)
mapping_index="path/to/minimap/index"   # path to index for mapping
tax_file="path/to/metadata"             # file containing taxonomy info
results_dir="metamap"                   # user-defined output directory for mapping results
strictness="--sensitive"                # options are --sensitive | --specific (in development)
experiment_id=""                        # will be used to name interactive output table

# Check for tools
if ! command -v samtools &> /dev/null; then
    echo "ERROR: samtools is not installed. Please install it and try again."
    exit 1
elif ! command -v minimap2 &> /dev/null; then
    echo "ERROR: minimap2 is not installed. Please install it and try again."
    exit 1
elif ! command -v fastp &> /dev/null; then
    echo "ERROR: fastp is not installed. Please install it and try again."
    exit 1
elif ! command -v multiqc &> /dev/null; then
    echo "ERROR: multiqc is not installed. Please install it and try again."
    exit 1
elif ! command -v msamtools &> /dev/null; then
    echo "ERROR: msamtools is not installed. Please install it and try again."
    exit 1
elif ! command -v salmon &> /dev/null; then
    echo "ERROR: salmon is not installed. Please install it and try again."
    exit 1
elif ! command -v bowtie2 &> /dev/null; then
    echo "ERROR: msamtools is not installed. Please install it and try again."
    exit 1
elif ! command -v fastqc &> /dev/null; then
    echo "ERROR: fastqc is not installed. Please install it and try again."
    exit 1
fi

# Check options
if [[ "$strictness" != "--sensitive" ]] && [[ "$strictness" != "--specific" ]]; then
    echo "MAIN: ERROR: Invalid strictness option."
    echo "  - Options are --sensitive | --specific"
    echo "MAIN: Exiting"
    exit 1
fi

# Create log file
log_filename="log-$(basename $0 .sh)_$(date +"%Y%m%d_%H%M%S").log"
exec > >(tee -a "$log_filename") 2>&1

# Log settings
echo "----- Starting control script for $version -----"
echo "MAIN: Start time: $(date +"%Y"."%m"."%d %H":"%M":"%S")"
echo "MAIN: Run Settings"
echo "  - working_dir:         $(pwd)"
echo "  - starting_format:     ${starting_format}"
echo "  - starting_dir:        ${starting_dir}"
echo "  - contaminant_dir:     ${contaminant_dir}"
echo "  - mapping_index:       ${mapping_index}"
echo "  - tax_file:            ${tax_file}"
echo "  - results_dir:         ${results_dir}"
echo "  - strictness:          ${strictness}"
echo "   "

# QC
echo "MAIN: Staring QC script..."
bash subscripts/qc.sh "$starting_format" "$starting_dir" > /dev/null
echo "  "

# Contaminant filtering (human and univec)
echo "MAIN: Starting contaminant filtering script..."
bash subscripts/filter_contaminants.sh trimmed_fastq "$contaminant_dir" > /dev/null
echo "   "

# Map to reference database
echo "MAIN: Starting mapping script..."
bash subscripts/meta_map.sh host-filt_fastq "$mapping_index" initial_mapping_results "$strictness" > /dev/null
echo "   "

# Generate mapping stats
echo "MAIN: Starting stat summary script..."
bash subscripts/mapping_stat_summary.sh initial_mapping_results
echo "   "

# Filter accessions and subset fasta 
echo "MAIN: Starting filter and subset script..."
mkdir -p subset_fastas
Rscript subscripts/filter_and_subset_fasta.R initial_mapping_results/mapping_stats "$mapping_reference" "$strictness" "$tax_file"

# Consolidate subset fastas if script was successful
if [[ -z "$(find "initial_mapping_results/mapping_stats" -type f -name *_subset.fasta)" ]]; then
    echo "MAIN: WARNING: No subset fastas were found in initial_mapping_results/mapping_stats. Script error or consolidation was already performed."
else
    echo "MAIN: Found subset fastas. Consolidating in subset_fasta directory..."
    mv initial_mapping_results/mapping_stats/*subset.fasta subset_fastas
fi
echo "   "

# Remap reads to filtered subset
echo "MAIN: Remapping reads to filtered fastas..."
bash subscripts/map_filtered_reads.sh "$strictness" > /dev/null
echo "  "

# Run Salmon quantification on fasta subsets
echo "MAIN: Starting salmon quantification..."
bash subscripts/salmon_quantification.sh "filtered_mapping_results" "subset_fastas" > /dev/null
echo "   "

# Compile final stats
echo "MAIN: Compiling final stat summary..."
cp salmon_quants/all_sample_quants/* ./final_stats/
Rscript subscripts/final_stat_summary.R final_stats "$tax_file"
echo "   "

# Create report
echo "MAIN: Creating interactive report..."
Rscript subscripts/helper_scripts/create_report_table.R final_stats "$experiment_id"
echo "   "

echo "MAIN: All scripts run successfully."

echo "  - Consolidating outputs in "$results_dir"..."
mkdir -p "$results_dir"
mkdir -p "$results_dir/logs"

mv trimmed_fastq "$results_dir"
mv host-filt_fastq "$results_dir"
mv initial_mapping_results "$results_dir"
mv subset_fastas "$results_dir"
mv log-* "$results_dir"/logs
mv fastp_reports "$results_dir"
mv filtered_mapping_results "$results_dir"
mv salmon_quants "$results_dir"
mv final_stats "$results_dir"

rm "$results_dir"/final_stats/*.sf
rm -r *temp*
rm -r subscripts

if [[ -d fastq ]]; then
    rm -r fastq
fi

# Compress intermediates
echo "  - Compressing intermediate files and folders..."
tar -czf "$results_dir"/metamap_intermediates.tar.gz \
    "$results_dir"/trimmed_fastq \
    "$results_dir"/host-filt_fastq \
    "$results_dir"/subset_fastas \
    "$results_dir"/salmon_quants

# Remove uncompressed files
if [[ -s "$results_dir"/metamap_intermediates.tar.gz ]]; then
    echo "  - Removing uncompressed intermediates..."
    rm -r \
        "$results_dir"/trimmed_fastq \
        "$results_dir"/host-filt_fastq \
        "$results_dir"/subset_fastas \
        "$results_dir"/salmon_quants
fi   
echo "MAIN: End time: $(date +"%Y"."%m"."%d %H":"%M":"%S")"