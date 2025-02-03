# Final stat summary
suppressMessages(suppressWarnings(library("tidyverse")))

# Writing mapping summmary in R
args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments is provided
if (length(args) != 2) {
  stop("R: This script requires 2 arguments: \n  
  - $1: path to final stats folder with salmon quants and mapping summary \n  
  - $2: file with taxonomy info \n")
}

directory <- args[1]
tax_file <- args[2]

taxonomy <- read.delim(paste0(tax_file), sep = "\t") %>% rename(Accession = accession)

# check tax file
if (nrow(taxonomy) == 0) {
  cat("R: WARNING: taxonomy file is empty, will not link taxa to stats. \n")
}

cat("R: Starting final summary script. \n")
cat("  - Using stat file directory:", paste0(directory), "\n")

quant_files <- list.files(path = directory, pattern = "quant.sf", full.names = T)
quant_in <- quant_files[file.size(quant_files) > 0]

mapping_stats <- read.csv(paste0(directory, "/final_mapping_stats.csv"), header = T, strip.white = T)

# Check inputs
if (length(quant_in) == 0) {
  stop("R: ERROR: No salmon quantification files were found. Exiting.")
}
  
if (nrow(mapping_stats) == 0) {
  stop("R: ERROR: No data in mapping stats file. Exiting.")
}

# Read in quant data
cat("R: Reading in quant files... \n")
quants <- lapply(quant_in, function(file) {
  cat("  - Processing: ", file, "\n")
  data <- read.delim(file, header = TRUE)
  data$File_Origin <- file
  return(data)
}) %>%
  bind_rows() %>%
  mutate(
    Sample = str_remove(File_Origin, paste0(directory, "/")),
    Sample = str_remove(Sample, "_quant.sf")
    ) %>%
  select(-File_Origin, -Length, "Accession" = Name, "Salmon_TPM" = TPM, "Salmon_read_number" = NumReads)

# Add mapping stats
final_stats <- mapping_stats %>% 
  left_join(quants, by = c("Sample", "Accession")) %>%
  group_by(Sample) %>%
  filter(Reads_mapped >= 20) %>%
  filter(Salmon_TPM > 0) %>%
  filter(Per_cov_1x >= 800 | bases_covered_1x >= 40) %>%
  # Filter out cases of short, highly multimapped sequences
  filter(
    ifelse(bases_covered_1x <= 2000,
      percent_MM <= 80,
      percent_MM <= 98
    )) %>%
  # Filter out cases of low coverage, highly multimapped sequences
  filter(
    ifelse(Per_cov_1x <= 5,
      percent_MM <= 95,
      percent_MM <= 98
    ))
 
# Try to link taxonomy
if (nrow(taxonomy) > 0){
  cat("R: Found tax file.  Adding taxonomy paths... \n")
  final_stats <- final_stats %>% left_join(taxonomy, by = "Accession")
}

# Write out file
cat("R: Writing output file... \n")
write.table(
  final_stats,
  paste0(directory, "/metamap_final_stats.tsv"),
  col.names = T,
  row.names = F,
  quote = F,
  sep = "\t"
)

cat("R: Final stat summary completed successfully.  Summary can be found in", directory, "/metamap_final_stats.csv", "\n")