suppressMessages(suppressWarnings(library("tidyverse")))
options(warn = -1)

args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments is provided
if (length(args) != 1) {
  stop("ERROR: R: This script requires one argument. The directory containing all stat files.")
}

# Set path to files
directory <- args[1]
cat("Directory:", directory, "\n")

# Get files
cat("R: Reading in files", "\n")
cat("  - Reading in idx file", "\n")
idx_files <- list.files(path=paste(directory), pattern="idxstat.tsv", full.names = T)
idx_in <- idx_files[file.size(idx_files) > 0]

cat("  - Reading in coverage files", "\n")
cov_files <- list.files(path=paste(directory), pattern="coverage.tsv", full.names = T)
cov_in <- cov_files[file.size(cov_files) > 0]

cat("  - Reading in library sizes", "\n")
libs <- read.delim(paste0(directory, "/reads_per_library.txt"), header = F) %>%
  select("Sample" = 1, "Library_Size" = 2)

empty_files <- c(
  "  R: EMPTY idx files: ", gsub(paste0(directory, "/"),"", idx_files[file.size(idx_files) == 0]), "\n",
  "  R: EMPTY cov files: ", gsub(paste0(directory, "/"),"", cov_files[file.size(cov_files) == 0])
)

cat("R: The following files had no data to read in.", "\n", empty_files, "\n")

if (all(lengths(list(idx_in, cov_in)) == 0)) {
  stop("ERROR: R: No data found for mapping statistics, indicates no reads were successfully mapped to the reference for samples provided. Exiting.")
}

# Read in files and format
## Index stats
cat("R: Reading in idx stat files \n")
idx <- lapply(idx_in, function(file) {
  cat("  - Processing: ", file, "\n")
  data <- read.delim(file, header = TRUE)
  # Check that file is bigger than just a header
  if (nrow(data) > 0) {
    data$File_Origin <- file
    return(data)
  }
}) %>% bind_rows() %>%
  mutate(
    Sample = str_remove(Sample, paste0(directory)),
    Sample = gsub(paste0("_idxstat.tsv"), "", Sample),
    Sample = str_remove(Sample, "/")
  )

## Coverage stats
cat("R: Reading in coverage data (this may take some time) \n")
cov <- lapply(cov_in, function(file) {
  cat("  - Processing: ", file, "\n")
  data <- read.delim(file, header = FALSE)
  data$File_Origin <- file
  data <- subset(data, V3 > 0)
  return(data)
}) %>% bind_rows() %>%
  select("Accession" = 1, "Position" = 2, "Depth" = 3, "Sample" = 4) %>%
  mutate(
    Sample = str_remove(Sample, paste0(directory)),
    Sample = gsub(paste0("_coverage.tsv"), "", Sample),
    Sample = str_remove(Sample, "/"),
  ) %>%
  left_join(libs, by = "Sample") %>%
  group_by(Sample, Accession) %>%
  left_join(idx, by = c("Sample", "Accession")) %>%
  summarize(
    bases_covered_1x = round(sum(Depth >= 1)),
    bases_covered_10x = round(sum(Depth >= 10)),
    Per_cov_1x = round(sum(Depth >= 1) / Genome_Length * 100, digits = 2),
    Per_cov_10x = round(sum(Depth >= 10) / Genome_Length * 100, digits = 2),
    Reads_mapped = Reads_Mapped,
    Library_Size = Library_Size,
    Num_MM_reads = Num_MM_reads,
    percent_MM = percent_MM,
    unique_reads = unique_reads,
    percent_mapped = (Reads_mapped / Library_Size) * 100,
    Genome_Length = Genome_Length,
    mean_depth_mapped = round(mean(Depth), digits = 3),
    max_depth = max(Depth),
    median_depth = median(Depth),
    .groups = 'drop'
  ) %>%
  ungroup() %>%
  group_by(Sample) %>%
  distinct(Accession, .keep_all = T)

## Combine stats and write out
cat("R: Combining stats and writing output file 'stat_summary.csv' \n")

stats_combined <- cov %>%
  group_by(Sample) %>%
  mutate(
    RPK = (Reads_mapped / (Genome_Length / 1000)),
    RPK_Factor = sum(RPK) / 1000000,
    TPM = RPK/RPK_Factor) %>%
  select(-RPK, -RPK_Factor) %>%
  mutate(Accession = str_remove(Accession, ";$"))

write.csv(stats_combined, paste0(directory, "/stat_summary.csv"), row.names = F, quote = F)

cat("R: RScript finished successfully.  Stat summary file found in ", directory, "\n")