#!/bin/bash

set -eu

echo "Running $0"
echo 'Start time:' `date +'%Y-%m-%d %T'`
echo "  "

# Set variables
fastq=""
database="/path/to/pluspfp_database"
output="kraken2_results/plus_pfp"

# Log parameters
echo "----- Run settings -----"
echo "fastq = ${fastq}"
echo "database = ${database}"
echo "output = ${output}"
echo "   "

mkdir -p ${output}
mkdir -p ${output}/reports
mkdir -p ${output}/krona_plots
mkdir -p ${output}/read_assignments

echo "----- Starting kraken2 -----"
echo "Checking input files..."
for i in ${fastq}/*_R1.fastq.gz; do
    samp=$(basename $i _R1.fastq.gz)

    # Check input
    if [ -s $i ] && [ -s ${fastq}/${samp}_R2.fastq.gz ]; then
        echo "${samp} ready for classifying"
    else
        echo "ERROR: Issue with ${samp}. File does not exist or is empty. Exiting."
        echo "  R1: $i"
        echo "  R2: ${samp}_R2.fastq.gz"

        exit 1
    fi

    # Run kraken2
    echo "Classifying ${samp}"
    kraken2 --threads 10 \
        --db "$database" \
        --report ${output}/reports/${samp}_kraken2_report.out \
        --paired  $i ${fastq}/${samp}_R2.fastq.gz | gzip -c >  ${output}/read_assignments/${samp}_kraken2_output.out.gz

    # Check output of Kraken2
    if [ -s ${output}/read_assignments/${samp}_kraken2_output.out.gz ]; then
        echo "${samp} successfully classified.  Creating krona plot."
    else
        echo "ERROR: issue classifying ${samp}."
    fi

    # Create krona plots
    echo "Creating krona plot for ${samp}"
    ktImportTaxonomy -tax /nobackup/lab_bergthaler/Tools/Krona/KronaTools/ \
        -m 3 \
        -t 5 \
        ${output}/reports/${samp}_kraken2_report.out \
        -o ${output}/krona_plots/${samp}_krona_plot.html

    # Check output
    if [ -e ${output}/krona_plots/${samp}_krona_plot.html ]; then
        echo "Krona plot generated for ${samp}"
        echo "   "
    else
        echo "ERROR: issue generating krona plot for ${samp}."
        echo "   "
    fi
done

echo "----- Kraken2 classification and krona plotting finished successfully -----"