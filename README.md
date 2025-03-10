# comparative_ww_virome_analysis_manuscript
All code used in the analysis and figures for the manuscript "Comparative wastewater virome analysis with different enrichment methods"

# Recreation of results
## Running the primary mapping workflow
1. Clone repo and download raw sequence files from European Nucleotide Archive.
```
# Clone repo
git clone https://github.com/thorntonmr/comparative_ww_virome_analysis_manuscript

# Get raw sequences
https://www.ebi.ac.uk/ena/browser/view/PRJEB86516
```

2. Create environment and install dependencies
```
# Create conda env (or mamba) from file
conda create -f environment.yml
conda activate metamap

# Install R dependencies
Rscript analysis_and_figures/install_r_dependencies.R
```

3. Get requisite files
- For the publication the contaminant set was created using the following
    - Human reference assembly: `GCF_000001405.26`
    - UniVec database: `https://ftp.ncbi.nlm.nih.gov/pub/UniVec/`
    - Concatenate the fastas 
    - Build the index with `bowtie2-build <contam_reference> <index_prefix>`

- Download mapping database: 
`wget https://zenodo.org/records/7876309/files/DB_v2.0.2.tar.gz`
    - For more information on mapping database see the publication by [Tisza et al., 2023](https://www.nature.com/articles/s41467-023-42064-1)
    - The resulting directory will contain the `mapping_reference` fasta, `mapping_index`, and `tax_file` needed subsequently

4. Setup config block
```md
# open the main.sh script in a text editor (scripts/main.sh)

# Set Variables
starting_format="bam"                                       # input file format (muts be bam or fastq)
starting_dir="path/to/starting/dir"                         # directory path for raw input files
contaminant_dir="path/to/bt2_ref"                           # path to bowtie2 index for contaminant reference
mapping_reference="path/to/virus_pathogen_database.fna"     # path to reference sequence(s)
mapping_index="path/to/virus_pathogen_database.mmi"         # path to index for mapping
tax_file="path/to/virus_pathogen_database.all_metadata.tsv" # file containing taxonomy info
results_dir="metamap"                                       # user-defined output directory for mapping results
strictness="--sensitive"                                    # do not change
experiment_id=""                                            # will be used to name interactive output table
```

5. Run script 
`bash main.sh` - creates a directory with all output files. 
- Subsequent analyses for the NGS data will work with the `$results_dir/final_stats/metamap_final_stats.tsv` file
- Running this workflow will produce the count data used in the sequencing portion of the results. Specifically it will yield a file that contains read counts per taxa (`Salmon_read_number`) that forms the basis for most of the figure recreations. 

6. Run Kraken2 for superkingdom level sample composition
    - Download Kraken2 [PlusPFP database](https://benlangmead.github.io/aws-indexes/k2) (warning: large file requiring ~143 GB space)
    - Install [KronaTools](https://github.com/marbl/Krona/wiki/KronaTools) 
        - Set binary path `export PATH=$PATH:/path/to/Krona/KronaTools/bin`
    - Open `scripts/kraken2.sh` in text editor
        - Add path to fastq files: in the publication, the post-QC fastq files `<filename>_trimmed.fastq.gz` were used 
        - Add path to PlusPFP database  
    - Run `scripts/kraken2.sh`

7. Run `scripts/lofreq.sh` to produce `.vcf` files used for nucleotide diversity calculations.
    - As before open the file in a text editor and edit the appropriate fields in the `# Set Variables` block

## Recreating the figures
- All files and scripts needed to recreate the figures from the publication are contained in `analysis_and_figures` directory. Knitting the Rmarkdown `figure_draft.Rmd` will reproduce figures used in the publication from pre-processed data.
- To independently reproduce the figures from the output of the pipeline above, code would need to be adapted to reflect the appropriate file paths.
