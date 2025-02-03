# Create table and visualize report

# Check libs
required_packages <- c("reactable", "reactablefmtr", "viridis")
for (package_name in required_packages) {
  if (!require(package_name, character.only = TRUE)) {
    cat("R: WARNING: The package", package_name, "is not installed. Installing now...")
    install.packages(package_name, dependencies = TRUE)
  }
}

suppressMessages(suppressWarnings(library("tidyverse")))
suppressMessages(library("reactable"))
suppressMessages(library("reactablefmtr"))
suppressMessages(library("viridis"))

args <- commandArgs(trailingOnly = TRUE)

directory <- args[1]
experiment <- args[2]

cat("R: Creating report...", "\n")
cat("  - Reading in stat file",  "\n")

stats <- read.delim(paste0(directory, "/metamap_final_stats.tsv"), sep = "\t") %>%
  mutate(
    Sample = str_remove(Sample, "_trimmed_host_filtered_mapped"),
    #Sample = word(Sample, 1, -2, sep = "_"),
    percent_MM = round(percent_MM, digits = 2),
  ) %>%
  select(
    Sample,
    "Species" = species, 
    Accession,
    Salmon_TPM,
    Per_cov_1x,
    Per_cov_10x,
    percent_MM,
    unique_reads,
    bases_covered_1x,
    bases_covered_10x,
    Reads_mapped,
    max_depth,
    median_depth,
    Genome_Length,
    genus,
    family
  )

cat("  - Creating stat table",  "\n")
final_table <- stats %>%
  reactable(
    .,
    pagination = TRUE,
    filterable = TRUE,
    showPageSizeOptions = TRUE,
    pageSizeOptions = c(10, 20, 100),
    defaultPageSize = 10,
    columns = list(
      Per_cov_1x = colDef(
        cell = data_bars(
          data = .,
          fill_color = viridis::inferno(5, direction = -1),
          background = '#F1F1F1',
          min_value = 0,
          max_value = 100,
          round_edges = TRUE,
          text_position = 'outside-end'
        )
      ),
      Per_cov_10x = colDef(
        cell = data_bars(
          data = .,
          fill_color = viridis::inferno(5, direction = -1),
          background = '#F1F1F1',
          min_value = 0,
          max_value = 100,
          round_edges = TRUE,
          text_position = 'outside-end'
        )
      ),
      percent_MM = colDef(
        cell = data_bars(
          data = .,
          fill_color = viridis(5, direction = -1),
          background = '#F1F1F1',
          min_value = 0,
          max_value = 100,
          round_edges = TRUE,
          text_position = 'outside-end'
        )
      ),
      Salmon_TPM = colDef(
        maxWidth = 100, 
        style = color_scales(., colors = c("grey", "gold", "firebrick2"), bias = 2), 
        format = colFormat(digits = 4)
      )
    )) %>% 
  add_title(paste0(experiment, ": Virus Detection Overview")) %>%
  google_font(font_family = "Inter")

cat("  - Saving table...",  "\n")

save_reactable_test(
  final_table,
  paste0(directory, "/", experiment, "_stat_table.html")
  )