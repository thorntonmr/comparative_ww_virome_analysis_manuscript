# Wastewater metagenomics
The repository contains all analysis scripts for the study on wastewater concentration methods for viral sequencing.

## Dependencies
- fastqc
- multiqc
- samtools
- minimap2
- msamtools
- fastp
- salmon
- bowtie2

- R packages:
    - tidyverse
    - viridis
    - Biostrings
    - reactable
    - reactablefmtr
    - Rsamtools

Installation should work easily enough with conda

```
mamba create -n metamap -c bioconda samtools msamtools minimap2 fastp multiqc salmon bowtie2 fastqc
```

## Usage
Currently the easiest way to run the analysis is to copy the scripts to the experiment directory.
- `main.sh` is the control script.  The only user input needed its to edit the top config block to include your paths and folder names and everything will run from there.  This can simply be copied from this repo, all dependencies will automatically be copied to the directory with the main script.
    - **NOTE**: it is best to run `main.sh` in the directory where you would like to have the results.

### Config block

```
# Set Variables
starting_format="fastq"
starting_dir="/home/mthomt69/personal_folder_BACKUP/gitlab/wastewater-metagenomics/sim_data/sim_v1_mutations/sim_fastq"
contaminant_dir="/mnt/winshare/Projects/Resources/Reference_Seqs/univec_human" 
mapping_reference="/home/mthomt69/personal_folder_BACKUP/gitlab/databases/twist_DB_v2.0.2/virus_pathogen_database.fna"
mapping_index="/home/mthomt69/personal_folder_BACKUP/gitlab/databases/twist_DB_v2.0.2/virus_pathogen_database.mmi"
tax_file="/home/mthomt69/personal_folder_BACKUP/gitlab/databases/twist_DB_v2.0.2/virus_pathogen_database.all_metadata.tsv"
results_dir="metamap_sim"
strictness="--sensitive"
experiment_id="sim_data"
```

- `starting_format`: two options are accepted.  `bam` or `fastq`.  `bam` currently assumes input with CeMM BSF sample naming format (will look for the "#" to parse the sample name).  `fastq` option requires the files to be in *R1/2.fastq.gz format.
- `starting_dir`: Directory path to the folder with the .bam or .fastq.gz files.  Best practice to include absolute file paths here.
- `contaminant_dir`: path to bowtie2 index of contaminants reference.  Default is set to a collection of UniVec DB references and the human reference genome.
- `mapping_reference`: path to fasta file for reference genome.  Will not change if using EsViritu DB
- `mapping_index`: path to minimap2 `.mmi` index of reference genome.  Will not change if using EsViritu DB.
- `tax_file`: path to file containing taxonomy linked to accession number.  Currently expects a `.tsv` with the first column title "accession".  This usually will not change if running with the EsViritu database
- `results_dir`: User defined name of final output directory for all files/folders produced in the analysis
- `strictness`: default `--sensitive` though `--specific` option is in development. Sets parameters for filtering reads and reference genomes.
    - `--sensitive`: ANI_filter=90%, min_reads=20, min_cov_1x=800bp|40%, %_MM=98%
- `experiment_id`: will be used to title the interactive stats table

### Output
- All outputs will be directed to the specified `$results_dir` at the end of the script
    - metamap_intermediates.tar.gz contains compressed intermediate files/folders.
        - trimmed_fastq
        - host-filt_fastq
        - subset_fastas
        - salmon_quants
    - individual files or folders can be extracted if needed but relative path must be specified
    - Example:
        ```tar -xcvf metamap_intermediates.tar.gz ./$results_dir/host-filt_fastq```
- `./fastp_reports`: contains all qc reports from initial read qc -> multiqc_report.html is the summary report
- `./intial_mapping_stats`: contains stats and both filtered/unfiltered bams from the first "filtering" mapping step
    - stat_summary.csv: all sequence mapping stats
    - excluded_accessions.tsv: lists all excluded accessions and the reason(s) for exclusion
    - excluded_species_mm.tsv: lists all species where there were >0 accession(s) excluded due to multimappers only 
- `./filtered_mapping_results`: second round stats from re-mapping to filtered reference set
- `./final_stats`
    - final_mapping_stats.csv: mapping stats prior to second filtering 
    - metamap_final_stats.tsv = final output file (used to construct interactive table)
    - interactive table

## Versions
- v1.0
- v1.1:
    - Added conditional filtering steps for multimapping thresholds based on bases covered
    - if bases_covered_1x < 2500 then the MM threshold was set to 50%
- v1.2:
    - Change multimapper definition to MAPQ == 0 (previously used bitflag + differing qname + rname to call MM)
    - Added one more conditional filter to MM thresholds
    - Changed initial QC to increase the read quality stringency to >=25 (Phred score)
- v1.3:
    - Added fuzzy multimappers to filter criteria.  Fuzzy multimappers are reads which multimap to another reference but have an edit distance of 1 (ie. 1 mismatch) with respect to the best mapping reference site.  In this update these reads will all be considered equal MM by setting MAPQ == 0. 
    - Filtered reads to keep only pairs which map to the same reference.
    - Changed filtering of short sequences due to MM% from <50% to <80% for all seqs <2500bp & ref_length >=5000bp (based on spike-in influenza data)
    - Removed `--best-hit` option from first round read filtering. 
    - Added a second filtering round prior to reporting detections. 
        - Min reads >= 20
        - Salmon_TPM > 0 (simulated data showed this to increase precision)
        - Sequences with <800bp | 40% covered 1x removed
        - If covered length is <2000bp, the multimap limit is set to 80%
        - If the covered percent is <5%, the multimap limit is set to 95% 
    - Added output file with excluded accessions + reason

## Limitations
- This script is designed to emphasize sensitivity over specificity.  This comes at the expense of including false positive calls but should maximize the reduction in false negatives.  In simulated analyses most of the false positive calls are a result of highly similar segmented virus calls (namely norovirus and sapporovirus) which can be dealt with by collapsing stats on the species level.  However, it is still necessary to be critical of the results.
