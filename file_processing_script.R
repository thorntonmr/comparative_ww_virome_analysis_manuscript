# Data file processing for dPCR and sequencing data for WWDv_001-3
## dPCR files were exported from QIAcuity software and manually compiled in excel to form "dPCR_all_results_updated.csv"

setwd("analysis_and_figures")

process_vcf <- function(vcf_file) {
  input <- readLines(vcf_file)
  
  # Filter out lines starting with '##'
  vcf_trim <- input[!grepl("^##", input)]
  
  # write tab-sep df
  vcf_df <- read.table(text = vcf_trim, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  
  # name columns
  vcf_df <- vcf_df %>% 
    dplyr::rename(
      accession = V1, 
      ref_position = V2,
      id = V3,
      ref_allele = V4,
      alt_allele = V5,
      quality = V6,
      filter = V7,
      info = V8
    )
  
  # add file name for ID
  vcf_df$file_origin <- str_remove(vcf_file, paste0(vcf_dir, "\\/"))
  
  vcf_df <- vcf_df %>%
    separate(info, into = c("depth", "allele_frequency", "strand_bias", "counts"), sep = ";") %>% 
    mutate(
      depth = as.numeric(str_remove(depth, "DP=")),
      allele_frequency = as.numeric(str_remove(allele_frequency, "AF="))
    ) %>% 
    separate(counts, into = c("type_count", "values"), sep = "=") %>%
    separate(values, into = c("ref_fw", "ref_rv", "alt_fw", "alt_rv"), sep = ",") %>%
    rowwise() %>%
    # combine fw and rv counts to single measure
    mutate(
      ref_count = sum(as.numeric(ref_fw), as.numeric(ref_rv)), 
      alt_count = sum(as.numeric(alt_fw), as.numeric(alt_rv))
    ) %>%
    select(accession, ref_position, ref_allele, alt_allele, ref_count, alt_count, depth, file_origin)
  
  # Return the trimmed data
  return(vcf_df)
}
get_nt_div <- function(vcf_file) {
  input <- readLines(vcf_file)
  
  # Funciton to calculate pairwise nucleotide diversity
  get_di_num <- function(x){
    if (length(x) > 1) {
      product_vec <- combn(x, m = 2, prod)
      di <- sum(product_vec)
      return(di)
    } else{
      # If only 1 nucleotide is present, di is 0
      return(0)
    }
  } 
  
  cat("Processing file", vcf_file, "...\n")
  
  # Filter out lines starting with '##'
  vcf_trim <- input[!grepl("^##", input)]
  
  # write tab-sep df
  vcf_df <- read.table(text = vcf_trim, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  
  # name columns
  vcf_df <- vcf_df %>% 
    dplyr::rename(
      accession = V1, 
      ref_position = V2,
      id = V3,
      ref_allele = V4,
      alt_allele = V5,
      quality = V6,
      filter = V7,
      info = V8
    )
  
  # add file name for ID
  vcf_df$file_origin <- vcf_file
  
  cat("  - Counting nucleotides per position\n")
  vcf_df <- vcf_df %>%
    separate(info, into = c("depth", "allele_frequency", "strand_bias", "counts"), sep = ";") %>% 
    mutate(
      depth = as.numeric(str_remove(depth, "DP=")),
      allele_frequency = as.numeric(str_remove(allele_frequency, "AF="))
    ) %>% 
    separate(counts, into = c("type_count", "values"), sep = "=") %>%
    separate(values, into = c("ref_fw", "ref_rv", "alt_fw", "alt_rv"), sep = ",") %>%
    rowwise() %>%
    # combine fw and rv counts to single measure
    mutate(
      ref_count = sum(as.numeric(ref_fw), as.numeric(ref_rv)), 
      alt_count = sum(as.numeric(alt_fw), as.numeric(alt_rv))
    ) %>%
    select(accession, ref_position, ref_allele, alt_allele, ref_count, alt_count, depth, file_origin)
  
  
  # Pivot to long form with each accession having one row per position per nucleotide
  cat("  - Pivoting df to long format \n")
  vcf_df_long <- vcf_df %>%
    pivot_longer(cols = c("ref_allele", "alt_allele"), values_to = "nucleotide", names_to = "allele_type") %>% 
    pivot_longer(cols = c(alt_count, ref_count), names_to = "allele_count_type", values_to = "count") %>%
    # make sure non matching rows are removed (ex. where alt_allele_type = ref_count_type)
    mutate(
      allele_type = word(allele_type, 1, 1, sep = "_"),
      allele_count_type = word(allele_count_type, 1, 1, sep = "_"),
    ) %>%
    filter(allele_type == allele_count_type) %>%
    # Collapse duplicate reference entries (from when multiple alt alleles are present at a site)
    group_by(accession, ref_position, allele_type) %>%
    distinct(nucleotide, .keep_all = T) %>%
    select(accession, ref_position, nucleotide, count, depth, allele_type, file_origin) %>%
    filter(count > 0) # remove zero diversity sites (ie. where alt allele AF = 1)
  
  # Calculate nucleotide diverisity per site (Di)
  cat("  - Calculating pairwise nucleotide diversity \n")
  nuc_div <- vcf_df_long %>% 
    group_by(accession, ref_position, file_origin) %>% 
    reframe(
      di_num = get_di_num(count),
      di = di_num / ((depth^2 - depth) / 2)
    ) %>%
    select(-di_num) %>%
    group_by(accession, file_origin) %>%       
    distinct(ref_position, .keep_all = T) %>%
    select(accession, ref_position, di, file_origin)
  
  # Return the trimmed data
  return(nuc_div)
}

### Process dPCR and write output ----
dPCR_dat <- read.csv(paste0(working_dir, "data_files/dPCR_all_results_updated.csv")) %>%
  select("Sample" = 2, Target, "Concentration" = 7) %>%
  filter(
    Target != "HPV6", 
    !grepl("Direct_E|NSC|WSC|NTC", Sample)
  ) %>%
  mutate(
    Replicate = case_when(
      grepl("_E1", Sample) ~ "R1",
      grepl("_E2", Sample) ~ "R2",
      grepl("_E3", Sample) ~ "R3"
    ),
    Method = case_when(
      grepl("Liquid", Sample) ~ "VDC",
      grepl("Pellet", Sample) ~ "VDCP",
      grepl("PEG", Sample) ~ "PEG",
      grepl("Nano", Sample) ~ "NT",
      grepl("Ultra", Sample) ~ "UF",
      grepl("Membrane", Sample) ~ "MEM",
      grepl("Unconcentrated_E", Sample) ~ "UNC",
      grepl("_Spike", Sample) ~ "UNC_Sp",
      grepl("Direct", Sample) ~ "DE",
      .default = Sample
    ),
    Sample = gsub("VacMan_", "Vac. direct capture ", Sample),
    Sample = gsub(" Liquid", "", Sample),
    Sample = gsub("Direct_Extract", "Direct", Sample),
    Concentration_adj = case_when(
      Method == "NT" ~ Concentration * 4,
      grepl("UNC", Method) ~ Concentration * 80,
      TRUE ~ Concentration
    ),
    Target = ifelse(Target == "E1", "HPV-E1", Target),
    #Sample = ifelse(grepl("Vac", Sample), word(Sample, 1, 2, sep = "_"), word(Sample, 1, 1, sep = "_"))
    Sample = gsub("_E\\d", "", Sample)
  ) 

dPCR_dat <- read.csv(paste0(working_dir, "data_files/dPCR_all_results_updated.csv")) %>%
  select("Sample" = 2, Target, "Concentration" = 7) %>%    
  filter(grepl("Direct_Extract_1-10_Dilution", Sample)) %>% 
  filter(Target == "MS2" | Target == "LCMV") %>% 
  mutate(Concentration = 10*Concentration) %>% 
  mutate(
    Replicate = case_when(
      grepl("_E1", Sample) ~ "R1",
      grepl("_E2", Sample) ~ "R2",
      grepl("_E3", Sample) ~ "R3"
    ),
    Method = case_when(
      grepl("Direct", Sample) ~ "DE",
      .default = Sample
    ),
    Sample = gsub("Direct_1-10_Dilution", "Direct", Sample),
    Concentration_adj = case_when(
      Method == "NT" ~ Concentration * 4,
      grepl("UNC", Method) ~ Concentration * 80,
      TRUE ~ Concentration
    ),
    Target = ifelse(Target == "E1", "HPV-E1", Target),
    Sample = gsub("_E\\d", "", Sample)
  ) %>%
  rbind(dPCR_dat)

# Data in excel sheet already includes volume adujstment factors
tu_dpcr <- read.csv(paste0(working_dir, "data_files/sc2_quants_TU.csv")) %>%
  select(Sample, Target, 'Concentration' = concentration, "Concentration_adj" = concentration_adj) %>%
  mutate(
    Replicate = case_when(
      grepl("1", Sample) ~ "R1",
      grepl("2", Sample) ~ "R2",
      grepl("3", Sample) ~ "R3"
    ),
    Method = case_when(
      grepl("liquid", Sample) ~ "VDC",
      grepl("pellet", Sample) ~ "VDCP",
      grepl("PEG", Sample) ~ "PEG",
      grepl("Nano", Sample) ~ "NT",
      grepl("UF", Sample) ~ "UF",
      grepl("Mem", Sample) ~ "MEM",
      grepl("_spike", Sample) ~ "UNC_Sp",
      grepl("unconc.", Sample) ~ "UNC",
      grepl("Direct", Sample) ~ "DE",
      .default = Sample
    ),
    Sample = gsub("VacMan_", "Vac. direct capture ", Sample),
    Sample = gsub(" Liquid", "", Sample),
  ) %>%
  select(-Sample) %>%
  left_join(dPCR_dat %>% select(Sample, Method, Replicate) %>% group_by(Replicate) %>% distinct(Sample, .keep_all = T))

# Data in excel sheet already includes volume adujstment factors
tu_16S <- read.csv(paste0(working_dir, "data_files/16S_quants_TU.csv"), strip.white = T) %>%
  select(Sample, 'Concentration' = concentration_raw, "Concentration_adj" = concentration_adj_final) %>%
  group_by(Sample) %>%
  summarise(
    Concentration = mean(Concentration),
    Concentration_adj = mean(Concentration_adj) 
  ) %>%
  mutate(
    Target = "16S",
    Replicate = case_when(
      grepl("1", Sample) ~ "R1",
      grepl("2", Sample) ~ "R2",
      grepl("3", Sample) ~ "R3"
    ),
    Method = case_when(
      grepl("Liquid", Sample) ~ "VDC",
      grepl("Pellet", Sample) ~ "VDCP",
      grepl("PEG", Sample) ~ "PEG",
      grepl("Nano", Sample) ~ "NT",
      grepl("Ultra", Sample) ~ "UF",
      grepl("Membrane", Sample) ~ "MEM",
      grepl("_spike", Sample) ~ "UNC_Sp",
      grepl("Unconcentrated", Sample) ~ "UNC",
      grepl("Direct", Sample) ~ "DE",
      .default = Sample
    ),
    Sample = gsub("VacMan_", "Vac. direct capture ", Sample),
    Sample = gsub(" Liquid", "", Sample),
  ) %>%
  select(-Sample) %>%
  left_join(dPCR_dat %>% select(Sample, Method, Replicate) %>% group_by(Replicate) %>% distinct(Sample, .keep_all = T))

# Data in excel sheet already includes volume adujstment factors
tu_pmmov <- read.csv(paste0(working_dir, "data_files/pmmov_quants_TU.csv"), strip.white = T) %>%
  select(Sample, Target, "Concentration" = concentration, "Concentration_adj" = concentration_adj) %>%
  mutate(
    Replicate = case_when(
      grepl("1", Sample) ~ "R1",
      grepl("2", Sample) ~ "R2",
      grepl("3", Sample) ~ "R3"
    ),
    Method = case_when(
      grepl("liquid", Sample) ~ "VDC",
      grepl("pellet", Sample) ~ "VDCP",
      grepl("PEG", Sample) ~ "PEG",
      grepl("Nano", Sample) ~ "NT",
      grepl("UF", Sample) ~ "UF",
      grepl("Mem", Sample) ~ "MEM",
      grepl("_spike", Sample) ~ "UNC_Sp",
      grepl("unconc.", Sample) ~ "UNC",
      grepl("Direct", Sample) ~ "DE",
      .default = Sample
    ),
    Sample = gsub("VacMan_", "Vac. direct capture ", Sample),
    Sample = gsub(" Liquid", "", Sample),
    # Add factor 10 dilution from dPCR rerun of VDC and PEG samples by TU
    Concentration = ifelse(Method %in% c("VDC", "PEG"), Concentration * 10, Concentration)
  ) %>%
  select(-Sample) %>%
  filter(Method != "VDCP", !is.na(Replicate)) %>%
  left_join(dPCR_dat %>% select(Sample, Method, Replicate) %>% group_by(Replicate) %>% distinct(Sample, .keep_all = T))


# All together
dpcr <- rbind(dPCR_dat, tu_dpcr, tu_16S, tu_pmmov) %>%
  mutate(Viral_type = case_when(
    grepl("Vaccinia|Influenza|MS2|PhiX|LCMV", Target) ~ "Spike-in",
    TRUE ~ "Endogenous"
  ))

write.table(dpcr, "data_files/processed_dpcr_data.csv", sep = ",", quote = F, row.names = F)

### Process sequencing data and write output ----
seq_dat <- read.delim(paste0(working_dir, "data_files/metamap_final_stats.tsv"), sep = "\t") %>%
  rowwise() %>%
  mutate(
    Sample = str_remove(Sample, "_trimmed_host_filtered_mapped"),
    Sampleshort = unlist(str_split(Sample, "_"))[1],
    Experiment = ifelse(grepl("Hyb", Sample), "WWDv_002", "WWDv_003"),
    Replicate  = paste0("R",gsub("\\D", "", Sampleshort)),
    Replicate2 = case_when(
      grepl("1_", Sample) ~ "R1",
      grepl("2_", Sample) ~ "R2",
      grepl("3_", Sample) ~ "R3"
    ),
    WWTP = case_when(
      grepl("Vienna", Sample) ~ "WWTP I",
      grepl("Klost", Sample) ~ "WWTP II",
      grepl("Alland", Sample) ~ "WWTP III",
      TRUE ~ "WWTP I"
    ),
    Method = case_when(
      grepl("VML", Sample) ~ "VDC",
      grepl("PEG", Sample) ~ "PEG",
      grepl("NT", Sample) ~ "NT",
      grepl("UF", Sample) ~ "UF",
      grepl("MEM", Sample) ~ "MEM",
      grepl("UNC", Sample) ~ "UNC",
    ),
    Sample = gsub("VML", "VDC", Sample),
    Sample = word(Sample, 1, -2, sep = "_"),
    percent_MM = round(percent_MM, digits = 2)
  )

write.table(seq_dat, file = "data_files/processed_seq_data.tsv", sep = "\t", quote = F)



### Process kraken2 data and write output ----
kraken_dat = list()
kraken_in <- list.files("data_files/kraken2_results_plus_pfp", pattern = "report.out", full.names = T)

kraken_dat <- lapply(kraken_in, read_kraken) %>% 
  bind_rows %>% 
  mutate(Sample = word(Sample, 1, 1, sep = "_S1"))

domain_dat <- kraken_dat %>%
  filter(Tax_level %in% c("D", "U")) %>%
  mutate(
    Enrichment = ifelse(grepl("ctrl", Sample), "Control", "Hybrid-Capture"),
    Replicate = case_when(
      grepl("1", Sample) ~ "R1",
      grepl("2", Sample) ~ "R2",
      grepl("3", Sample) ~ "R3"
    ),
    Method = case_when(
      grepl("PEG", Sample) ~ "PEG",
      grepl("UF", Sample) ~ "UF",
      grepl("MEM", Sample) ~ "MEM",
      grepl("NT", Sample) ~ "NT",
      grepl("VML", Sample) ~ "VDC",
      grepl("UNC", Sample) ~ "UNC",
      grepl("DE", Sample) ~ "DE",
    ),
    Experiment = ifelse(grepl("Hyb", Sample), "WWDv_002", "WWDv_003"),
    WWTP = case_when(
      grepl("Alland", Sample) ~ "Alland",
      grepl("Kloster", Sample) ~ "Klosterneuburg",
      grepl("Vienna", Sample) ~ "Vienna",
      grepl("Hyb", Sample) ~ "Vienna"
    )
  ) 

write.table(kraken_dat, file = "data_files/processed_kraken_data.tsv", sep = "\t", quote = F)
write.table(domain_dat, file = "data_files/processed_domain_data.tsv", sep = "\t", quote = F)


### Process vcf data ----
vcf_dir <- "data_files/vcf_files"
vcf_list <- list.files(paste0(vcf_dir), pattern = "lofreq.vcf", full.names = T)
vcf <- lapply(vcf_list, get_nt_div) %>% bind_rows()

pi_dat <- vcf %>%
  mutate(
    Sample = word(file_origin, -1, -1, sep = "\\/"),
    Sample = str_remove(Sample, "_lofreq.vcf"),
    Experiment = ifelse(grepl("Hyb", Sample), "WWDv_002", "WWDv_003"),
    Replicate = case_when(
      grepl("1_", Sample) ~ "R1",
      grepl("2_", Sample) ~ "R2",
      grepl("3_", Sample) ~ "R3"
    ),
    WWTP = case_when(
      grepl("Vienna", Sample) ~ "WWTP 2",
      grepl("Klost", Sample) ~ "WWTP 3",
      grepl("Alland", Sample) ~ "WWTP 4",
      TRUE ~ "WWTP 1"
    ),
    Method = case_when(
      grepl("VML", Sample) ~ "VDC",
      grepl("PEG", Sample) ~ "PEG",
      grepl("NT", Sample) ~ "NT",
      grepl("UF", Sample) ~ "UF",
      grepl("MEM", Sample) ~ "MEM",
      grepl("UNC", Sample) ~ "UNC",
    ),
    Sample = gsub("VML", "VDC", Sample),
    Sample = word(Sample, 1, -2, sep = "_"),
  )

write.table(
  pi_dat, 
  "data_files/processed_individual_nt_diversity_file.tsv", 
  sep = "\t", 
  quote = F, 
  col.names = T, 
  row.names = F
  )
