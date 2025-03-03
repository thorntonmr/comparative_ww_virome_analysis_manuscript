# comparative_ww_virome_analysis_manuscript
All code used in the analysis and figures for the manuscript "Comparative wastewater virome analysis with different enrichment methods"

# Recreation of results
1. Clone repo and download raw sequence files from European Nucleotide Archive.
```
# Clone repo
git clone https://github.com/thorntonmr/comparative_ww_virome_analysis_manuscript

# Get sequences
```

2. Run sequence processing pipeline to get mapping results
- For this analysis we map to a reference set that can be downloaded from ``

2.1 Create environment
```
# Create conda env (or mamba) from file
conda create -f environment.yml
conda activate metamap
```

2.2 Get requisite file
- For the publication the contaminant set was created using the following
    - Human reference assembly: `GCF_000001405.26`
    - UniVec database: `https://ftp.ncbi.nlm.nih.gov/pub/UniVec/`

- Download mapping database: 
`wget https://zenodo.org/records/7876309/files/DB_v2.0.2.tar.gz`
    - For more information on mapping database see the publication by [Tisza et al., 2023](https://www.nature.com/articles/s41467-023-42064-1)
    - The resulting directory will contain the `mapping_reference` fasta, `mapping_index`, and `tax_file` needed subsequently

2.3 Setup config block
```md
# open the main.sh script in a text editor

# Set Variables
starting_format="bam"                                       # input file format (must be either .bam or _R*.fastq.gz)
starting_dir="path/to/starting/dir"                         # directory path for raw input files (bam or fastq)
contaminant_dir="path/to/bt2_ref"                           # path to bowtie2 index for contaminant reference
mapping_reference="path/to/virus_pathogen_database.fna"     # path to reference sequence(s)
mapping_index="path/to/virus_pathogen_database.mmi"         # path to index for mapping
tax_file="path/to/virus_pathogen_database.all_metadata.tsv" # file containing taxonomy info
results_dir="metamap"                                       # user-defined output directory for mapping results
strictness="--sensitive"                                    # options are --sensitive | --specific (in development)
experiment_id=""                                            # will be used to name interactive output table
```
