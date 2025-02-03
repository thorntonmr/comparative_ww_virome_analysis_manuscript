#!/bin/bash

set -eu

# Call SNPs in metagenomics data for downstream diversity assessment
# https://csb5.github.io/lofreq/

# Create log file
log_filename="lofreq.log"
exec > >(tee -a "$log_filename") 2>&1

# Set variables 
bam_dir="path/to/bams"
fasta_file="path/to/reference"
out_dir="vcf_files"
log="lofreq.log"

# Log variables
echo "working_dir:  $(pwd)"
echo "bam_dir:      $bam_dir" 
echo "fasta_file:   $fasta_file"
echo "out_dir:      $out_dir" 

mkdir -p $out_dir

# Run lofreq
echo "Running LoFreq..."
for i in "${bam_dir}"/*_mapped.bam; do
    samp=$(basename $i _trimmed_host_filtered_mapped.bam)
    
    echo "  - Processing ${samp}"
    lofreq call -f $fasta_file -o ${out_dir}/${samp}_lofreq.vcf "$i"
done

echo "Script finished successfully"