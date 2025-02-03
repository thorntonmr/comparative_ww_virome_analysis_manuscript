# Subsetting fasta for round 2 mapping
suppressMessages(suppressWarnings(library("tidyverse")))
suppressMessages(library(Biostrings))

# Writing mapping summmary in R
args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments is provided
if (length(args) != 4) {
  stop("R: This script requires two arguments: \n
    - $1 as the stat file path \n
    - $2 is the path to the sequence collection for subsetting \n
    - $3 is the stringency setting \n
    - $4 is the taxonomy file. \n"
    )
}

directory <- args[1]
full_fasta <- readDNAStringSet(args[2])
strictness <- args[3]
tax_file <- args[4]

cat("R: Starting script to subset fastas for preliminary detected viruses.", "\n")
cat("R: Using strictness setting:", strictness, "\n")

# Read in stats file
stat_file <- read.csv(paste0(directory, "/stat_summary.csv"), header = T, strip.white = T)
tax_table <- read.delim(tax_file, sep = "\t") %>% dplyr::rename(Accession = accession)

if (strictness == "--specific") {
  min_reads <- 10
  per_cov_filter <- 50
  base_cov_filter <- 1000
  per_mm_filter <- 98
  #reads_mm_filter <- 10

} else if (strictness == "--sensitive") {
  min_reads <- 20
  per_cov_filter <- 40
  base_cov_filter <- 800
  per_mm_filter <- 98
  #reads_mm_filter <- 20

} else {
  stop("Invalid strictness value provided")
}

cat("  - Preset options:", "\n")
cat("    - min reads:", min_reads, "\n")
cat("    - % cov filter:", per_cov_filter, "\n")
cat("    - base cov filter:", base_cov_filter, "\n")
cat("    - % multimapper filter:", per_mm_filter, "\n")
#cat("    - read number mm filter:", reads_mm_filter, "\n")

# Filter based on pre-determined criteria
filtered_list <- stat_file %>%
  group_by(Sample) %>%
  filter(Reads_mapped >= min_reads) %>%
  filter(Per_cov_1x >= per_cov_filter | bases_covered_1x >= base_cov_filter) %>%
  # Filter out cases of short, highly multimapped sequences
  filter(
    ifelse(Genome_Length >= 5000 & bases_covered_1x <= 2500,
      percent_MM <= 80,
      percent_MM <= per_mm_filter
    )) %>%
  # Filter out cases of short, highly multimapped sequences
  filter(
    ifelse(Per_cov_1x <= 5,
      percent_MM <= 95,
      percent_MM <= per_mm_filter
    )) %>%
  distinct(Accession)

# Get filtered species and add reason
excluded_acc <- stat_file %>%
  mutate(
    RN = ifelse(Reads_mapped < min_reads, "read_number", ""),
    LC = ifelse(bases_covered_1x < base_cov_filter & Per_cov_1x < per_cov_filter, "low_coverage", ""),
    MM = ifelse(percent_MM > per_mm_filter, "multimappers", ""),
    LCMM = ifelse(bases_covered_1x < 2500 & percent_MM > 80, "low_bp_cov_MM", ""),
    LGMM = ifelse(Per_cov_1x < 5 & percent_MM > 95, "low_per_cov_MM", "")
  ) %>%
  mutate(filter_reason = str_trim(paste(RN, LC, MM, LCMM, LGMM))) %>%
  mutate(filter_reason = ifelse(filter_reason == "", "pass", filter_reason)) %>%
  mutate(filter_reason = gsub(" ", ",", filter_reason)) %>%
  left_join(tax_table, by = "Accession")

summary_stats <- excluded_acc %>%
  group_by(Sample, species) %>%
  summarize(
    n_pass = sum(filter_reason == "pass"),
    n_failed_mm = sum(filter_reason == "multimappers"),
    n_failed_low_cov_mm = sum(filter_reason == "low_per_cov_MM"),
    .groups = 'drop'
  ) %>%
  filter(n_pass == 0 & (n_failed_mm > 0 | n_failed_low_cov_mm))

out_table <- excluded_acc %>%
  select(
    Sample,
    Accession,
    species,
    genus,
    family,
    bases_covered_1x,
    bases_covered_10x,
    Per_cov_1x,
    Per_cov_10x,
    percent_MM,
    Num_MM_reads,
    unique_reads,
    Reads_mapped,
    max_depth,
    filter_reason
  ) %>%
  filter(filter_reason != "pass")

# write out tables
cat("R: Writing records of excluded accessions...")
write.table(
  out_table,
  paste0(directory, "/excluded_accession_table.tsv"),
  col.names = T,
  row.names = F,
  quote = F,
  sep = "\t"
)

write.table(
  summary_stats,
  paste0(directory, "/excluded_mm_species_summary.tsv"),
  col.names = T,
  row.names = F,
  quote = F,
  sep = "\t"
)

# Subset fasta
subset_fasta_list <- list()
for (i in unique(filtered_list$Sample)) {
  sample <- filtered_list %>% filter(Sample == i) %>% pull(Accession)
  cat("\n")
  cat("R: Subsetting fasta for sample:",  i, "\n")
  cat("  - Found", length(sample), "accessions passing filter.", "\n")
  
  writeXStringSet(
    x = full_fasta[which(names(full_fasta) %in% sample)], 
    filepath = paste0(directory, "/", i, "_fasta_subset.fasta"), 
    format = "fasta"
  )
}

cat("\n")
cat("R: Subset script finished successfully.", "\n")