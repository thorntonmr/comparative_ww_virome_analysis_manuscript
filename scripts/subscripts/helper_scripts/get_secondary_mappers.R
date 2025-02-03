# Mapping stats
suppressMessages(suppressWarnings(library("tidyverse")))
suppressMessages(library("Rsamtools"))
options(warn = -1)

# Writing mapping summmary in R
args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments is provided
if (length(args) != 2) {
  stop("    ERROR: R: This script requires two arguments: $1 as the bam file path, $2 as the sample prefix")
}

bam_dir = args[1]
sample_prefix = args[2]
bam_file = paste0(bam_dir, "/", sample_prefix, ".bam")

cat("    R: Starting stat collection for", sample_prefix, "\n")
cat("      - Input options:", "\n")
cat("        - Bam file:", bam_file, "\n")

# Get idx stats
idx <- idxstatsBam(
  file = bam_file,
  index = paste0(bam_file, ".bai")
)

# Filter for secondary mapping reads
bam_flags_to_filter <- scanBamFlag(isSecondaryAlignment = TRUE)
info_to_extract <- c("qname", "flag", "rname", "pos", "mapq", "cigar", "mrnm", "mpos", "isize", "seq", "qual")

bam_data <- scanBam(
  file = bam_file,
  index = paste0(bam_file, ".bai"),
  param = ScanBamParam(
    #flag = bam_flags_to_filter,
    what = info_to_extract
      )
  )

# Convert bam view into dataframe
bam_info <- as.data.frame(bam_data[[1]])

# Summarize number of intergenomic multimappers per accession
mm_acc <- bam_info %>%
  filter(mapq == 0) %>%
  select(qname, "seqnames" = rname) %>%
  group_by(seqnames) %>%
  summarise(num_mm_reads = n()) %>%
  right_join(idx, by = "seqnames") %>%
  select("accession" = seqnames, seqlength, mapped, num_mm_reads, -unmapped) %>%
  mutate(
    num_mm_reads = ifelse(is.na(num_mm_reads), 0, num_mm_reads),
    percent_mm = num_mm_reads / mapped * 100,
    unique_reads = mapped - num_mm_reads, 
    Sample = sample_prefix
    ) %>%
  filter(mapped > 0) %>%
  dplyr::rename(
    Sample = Sample,
    Accession = accession,
    Genome_Length = seqlength,
    Reads_Mapped = mapped,
    Num_MM_reads = num_mm_reads,
    percent_MM = percent_mm,
    unique_reads = unique_reads
  )

# Write output file
cat("      - Writing file:", paste0(bam_dir, "/mapping_stats", sample_prefix, "_idxstat.tsv", "\n"))
write.table(
  mm_acc, 
  file = paste0(bam_dir, "/mapping_stats/", sample_prefix, "_idxstat.tsv"), 
  sep = "\t", 
  quote = F,
  col.names = T, 
  row.names = F
  )

cat("    R: Script finished successfully. Output files can be found in", bam_dir, "/mapping_stats", "\n")