# Input ----
library("tidyverse")
library("performance")
library("viridis")
library("lmtest")
library("patchwork")
library("marginaleffects")

conflicted::conflicts_prefer(dplyr::select)
conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::left_join)

auxdata <- read.delim("./analysis_and_figures/data_files/virus_genus_characteristic.tsv", sep = "\t") %>% 
  rowwise() %>%
  mutate(mean_size = round(mean(c(size_min, size_max)))) %>%
  mutate(sizeclass = case_when(
    mean_size >= 150 ~ "large", 
    dplyr::between(mean_size, 75, 150) ~ "mid", 
    mean_size <= 75 ~ "small")
  ) 

# Prep data ----
vir <- wwdv_002 %>%
  group_by(Method, Replicate, subspecies, species, genus, family, class) %>%
  summarize(Salmon_read_number = sum(Salmon_read_number)) %>%
  ungroup() %>%
  select(Method, Replicate, subspecies, species, genus, family, class, Salmon_read_number) %>%
  group_by(Method, subspecies) %>%
  mutate(mean_reads = mean(Salmon_read_number))

# Use mean of 3 unconcentrated samples as reference
ref <- vir %>% 
  filter(Method == "UNC") %>%
  mutate(ref_reads = mean_reads) %>%
  select(subspecies, species, genus, family, class, ref_reads) %>%
  distinct(subspecies, .keep_all = T)

# Divide each sample replicate by the mean of the reference
vir_joined <- vir %>% 
  dplyr::left_join(ref %>% ungroup() %>% select(-Method)) %>% 
  select(-mean_reads)

try <- vir_joined 

# define factor order
try$Method <- factor(try$Method, levels = c("UNC", "MEM", "UF", "NT", "VDC", "PEG"))
try <- try %>% left_join(auxdata)

# Model
lm2 <- lm(
  log10(Salmon_read_number / ref_reads) ~ Method + envelope + sizeclass + Nuc + Method:envelope + Method:sizeclass + Method:Nuc, 
  data = try
  )
summary(lm2)
performance::check_model(lm2)

# Plots for enrichment predictions morpho model
s1 <- plot_predictions(lm2, by = c("Method", "sizeclass")) + theme_bw() +  labs(title = "Estimates [log10]") + 
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.4)

s2 <- plot_predictions(lm2, by = c("Method", "envelope")) + theme_bw() +labs(title = "Estimates [log10]") + 
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.4)

s3 <- plot_predictions(lm2, by = c("sizeclass", "Method")) + theme_bw() + labs(title = "Estimates [log10]") + 
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.4)

s4 <- plot_predictions(lm2, by = c("envelope", "Method")) + theme_bw() +labs(title = "Estimates [log10]") + 
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.4)

(s1+s2) / (s3+s4) + plot_layout(guides = "collect")

# 1. Size comparisons ----
size_comparison <- comparisons(
  lm2,
  variables = "Method",
  by = "sizeclass"
) 

size_comparison <- size_comparison %>% mutate(contrast = str_remove_all(contrast, "mean\\(|\\)"))
size_comparison <- size_comparison %>% mutate(signif = case_when(
  p.value > 0.05 ~ "ns",
  p.value <= 0.001 ~ "sig***",
  p.value <= 0.01 ~ "sig**",
  p.value <= 0.05 ~ "sig*",
))

# Extract pairwise methods significance tests for each sizeclass level
out <- list()
for (i in 1:length(unique(try$sizeclass))) {
  size <- unique(try$sizeclass)[i]
  
  c <- comparisons(
    lm2, 
    variables = "Method", 
    newdata = datagrid(sizeclass = size), 
    hypothesis = "pairwise") %>% 
    as.data.frame() %>% 
    mutate(sizeclass = size) 
  
  out[[i]] <- c
}

out_size <- out %>% bind_rows() %>% group_by(sizeclass) %>% distinct(term, .keep_all = T)
adj.p <- p.adjust(out_size$p.value, method = "holm")

# create df to plot significance levels between methods within a sizeclass
max_val <- size_comparison %>% group_by(sizeclass) %>% summarize(max = max(conf.high) + 0.3)
sig <- out_size %>%
  select(term, sizeclass, p.value) %>%
  separate(term, into = c("group1", "group2"), sep = "\\) - \\(") %>%
  ungroup() %>%
  mutate(
    group1 = str_remove(group1, "\\("),
    group2 = str_remove(group2, "\\)"),
    p.adj = adj.p
  ) %>%
  mutate(sig = case_when(
    p.adj <= 0.05 & p.adj >= 0.01 ~ "*",
    p.adj < 0.01 & p.adj >= 0.001 ~ "**",
    p.adj < 0.001 ~ "***",
    p.adj > 0.05 ~ "ns"
  )) %>% 
  left_join(max_val) %>%
  filter(group1 == "PEG - UNC" | group2 == "PEG - UNC")

# Add raw fold enrichments per species per sizeclass as jitter
pdf <- try %>% 
  filter(Method != "UNC") %>%
  mutate(log10fc = log10(fold_enrichment_reads)) %>% 
  #left_join(count_fam) %>%
  mutate(contrast = paste0(Method, " - UNC"))

plot_c <- size_comparison%>% filter(term == "Method") %>% group_by(sizeclass) %>% distinct(contrast, .keep_all = T)

pdf$contrast <- factor(pdf$contrast, levels = c("PEG - UNC", "VDC - UNC", "UF - UNC", "NT - UNC", "MEM - UNC"))
plot_c$contrast <- factor(plot_c$contrast, levels = c("PEG - UNC", "VDC - UNC", "UF - UNC", "NT - UNC", "MEM - UNC"))

error_size <- pdf %>% 
  group_by(sizeclass, Method) %>% 
  summarize(q975 = quantile(log10fc, 0.975), q25 = quantile(log10fc, 0.025)) %>% 
  left_join(pdf %>% ungroup() %>% select(Method, contrast)) %>%
  distinct(contrast, .keep_all = T)

plot_c2 <- plot_c %>% left_join(error_size)

### Plot enrichment performance by each size level ----
plot_size <- ggplot(plot_c2, aes(x = contrast, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_violin(
    data = pdf %>% left_join(error_size) %>% group_by(Method, sizeclass) %>% filter(dplyr::between(log10fc, q25, q975)),
    aes(x = contrast, y = log10fc),
    fill = "grey",
    color = "grey75",
    alpha = 0.75
  ) +
  geom_jitter(
    data = pdf %>% left_join(error_size) %>% group_by(Method, sizeclass) %>% filter(!dplyr::between(log10fc, q25, q975)), 
    inherit.aes = F,
    aes(x = contrast, y = log10fc),
    alpha = 0.5,
    size = 0.5,
    width = 0.05
  ) +
  geom_errorbar(aes(ymin = q25, ymax = q975), width = 0.25, alpha = 0.75) +
  geom_point(shape = 21, alpha = 0.9, size = 3.5, aes(fill = signif)) +
  scale_fill_manual(values = c("grey", "#FEC368", "#FB8A30", "#F45C10")) +
  stat_pvalue_manual(
    data = sig %>% filter(sig != "ns"),
    label = "sig",
    xmin = "group1",
    xmax = "group2",
    y.position = "max",
    step.increase = 0.075,
    step.group.by = "sizeclass",
    coord.flip = F,
    bracket.size = 0.3,
    tip.length = 0.02,
    label.size = 3,
    bracket.nudge.y = 0.75
  ) +
  theme_bw() + 
  facet_wrap(~sizeclass, ncol = length(unique(plot_c$sizeclass)), strip.position = "top") + 
  labs(
    title = "Virion size",
    x = "",
    y = "Estimate [log10(fc)]",
    size = "Num. subspecies",
    fill = "Sig. level \n[method/unconc]"
  ) +
  theme(
    #plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8)
  )

# 2. Enveloped comparison ----
env_comparison <- comparisons(lm2, variables = "Method", by = "envelope") 

env_comparison <- env_comparison %>% 
  mutate(contrast = str_remove_all(contrast, "mean\\(|\\)")) %>% 
  mutate(signif = case_when(
    p.value > 0.05 ~ "ns",
    p.value <= 0.001 ~ "sig***",
    p.value <= 0.01 ~ "sig**",
    p.value <= 0.05 ~ "sig*"))

# Extract pairwise methods significance tests for each sizeclass level
out <- list()
for (i in 1:length(unique(try$envelope))) {
  class <- unique(try$envelope)[i]
  
  c <- comparisons(lm2, variables = "Method", newdata = datagrid(envelope = class), hypothesis = "pairwise") %>%
    as.data.frame() %>%
    mutate(envelope = class)
  
  out[[i]] <- c
}

out_env <- out %>% bind_rows() %>% group_by(envelope) %>% distinct(term, .keep_all = T)
adj.p <- p.adjust(out_env$p.value, method = "holm")

# create df to plot significance levels between methods within a sizeclass
max_val <- env_comparison %>% group_by(envelope) %>% summarize(max = max(conf.high) + 0.3)
sig <- out_env %>%
  select(term, envelope, p.value) %>%
  separate(term, into = c("group1", "group2"), sep = "\\) - \\(") %>%
  ungroup() %>%
  mutate(
    group1 = str_remove(group1, "\\("),
    group2 = str_remove(group2, "\\)"),
    p.adj = adj.p
  ) %>%
  mutate(sig = case_when(
    p.adj <= 0.05 & p.adj >= 0.01 ~ "*",
    p.adj < 0.01 & p.adj >= 0.001 ~ "**",
    p.adj < 0.001 ~ "***",
    p.adj > 0.05 ~ "ns"
  )) %>% 
  left_join(max_val) %>%
  filter(group1 == "PEG - UNC" | group2 == "PEG - UNC")

plot_e <- env_comparison %>% filter(term == "Method") %>% group_by(envelope) %>% distinct(contrast, .keep_all = T)

error_env <- pdf %>% 
  group_by(envelope, Method) %>% 
  summarize(q975 = quantile(log10fc, 0.975), q25 = quantile(log10fc, 0.025)) %>% 
  left_join(pdf %>% ungroup() %>% select(Method, contrast)) %>%
  distinct(contrast, .keep_all = T)

plot_e2 <- plot_e %>% left_join(error_env)

### Plot envelope effects ----
plot_env <- ggplot(plot_e2, aes(x = contrast, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_violin(
    data = pdf %>% left_join(error_env) %>% group_by(Method, envelope) %>% filter(dplyr::between(log10fc, q25, q975)),
    aes(x = contrast, y = log10fc),
    fill = "grey",
    color = "grey75",
    alpha = 0.75
  ) +
  geom_jitter(
    data = pdf %>% left_join(error_env) %>% group_by(Method, envelope) %>% filter(!dplyr::between(log10fc, q25, q975)), 
    inherit.aes = F,
    aes(x = contrast, y = log10fc),
    alpha = 0.5,
    size = 0.5,
    width = 0.05
  ) +
  geom_errorbar(aes(ymin = q25, ymax = q975), width = 0.25, alpha = 0.75) +
  geom_point(shape = 21, alpha = 0.9, size = 3.5, aes(fill = signif)) +
  scale_fill_manual(values = c("grey", "#FEC368", "#FB8A30", "#F45C10")) +
  stat_pvalue_manual(
    data = sig %>% filter(sig != "ns"),
    label = "sig",
    xmin = "group1",
    xmax = "group2",
    y.position = "max",
    step.increase = 0.075,
    step.group.by = "envelope",
    coord.flip = F,
    bracket.size = 0.3,
    tip.length = 0.02,
    label.size = 3,
    bracket.nudge.y = 0.75
  ) +
  theme_bw() + 
  facet_wrap(~envelope, ncol = length(unique(env_comparison$envelope)), strip.position = "top") + 
  labs(
    title = "Envelope status",
    x = "",
    y = "",
    size = "Num. subspecies",
    fill = "Sig. level \n[method/unconc]"  
  ) +
  theme(
    #plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
    legend.position = "none"
  )

# 3. Nucleic acid type ----
nuc_comp <- comparisons(
  lm2,
  variables = "Method",
  by = "Nuc"
) 

nuc_comp <- nuc_comp %>% 
  mutate(contrast = str_remove_all(contrast, "mean\\(|\\)")) %>% 
  mutate(signif = case_when(
    p.value > 0.05 ~ "ns",
    p.value <= 0.001 ~ "sig***",
    p.value <= 0.01 ~ "sig**",
    p.value <= 0.05 ~ "sig*"))

# Extract pairwise methods significance tests for each sizeclass level
out <- list()
for (i in 1:length(unique(try$Nuc))) {
  class <- unique(try$Nuc)[i]
  
  c <- comparisons(lm2, variables = "Method", newdata = datagrid(Nuc = class), hypothesis = "pairwise") %>%
    as.data.frame() %>%
    mutate(Nuc = class)
  
  out[[i]] <- c
}

out_nuc <- out %>% bind_rows() %>% group_by(Nuc) %>% distinct(term, .keep_all = T)
adj.p <- p.adjust(out_nuc$p.value, method = "holm")

# create df to plot significance levels between methods within a sizeclass
max_val <- nuc_comp %>% group_by(Nuc) %>% summarize(max = max(conf.high) + 0.3)
sig <- out_nuc %>%
  select(term, Nuc, p.value) %>%
  separate(term, into = c("group1", "group2"), sep = "\\) - \\(") %>%
  ungroup() %>%
  mutate(
    group1 = str_remove(group1, "\\("),
    group2 = str_remove(group2, "\\)"),
    p.adj = adj.p
  ) %>%
  mutate(sig = case_when(
    p.adj <= 0.05 & p.adj >= 0.01 ~ "*",
    p.adj < 0.01 & p.adj >= 0.001 ~ "**",
    p.adj < 0.001 ~ "***",
    p.adj > 0.05 ~ "ns"
  )) %>% 
  left_join(max_val) %>%
  filter(group1 == "PEG - UNC" | group2 == "PEG - UNC")

plot_n <- nuc_comp %>% filter(term == "Method") %>% group_by(Nuc) %>% distinct(contrast, .keep_all = T)

error_nuc <- pdf %>% 
  group_by(Nuc, Method) %>% 
  summarize(q975 = quantile(log10fc, 0.975), q25 = quantile(log10fc, 0.025)) %>% 
  left_join(pdf %>% ungroup() %>% select(Method, contrast)) %>%
  distinct(contrast, .keep_all = T)

plot_n2 <- plot_n %>% left_join(error_nuc)

### Plot nucleic acid effects ----
plot_nuc <- ggplot(plot_n2, aes(x = contrast, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_violin(
    data = pdf %>% left_join(error_nuc) %>% group_by(Method, Nuc) %>% filter(dplyr::between(log10fc, q25, q975)), 
    aes(x = contrast, y = log10fc),
    fill = "grey",
    color = "grey75",
    alpha = 0.75
    #flip = T,
    #position = position_nudge(x = -0.05)
  ) +
  geom_jitter(
    data = pdf %>% left_join(error_nuc) %>% group_by(Method, Nuc) %>% filter(!dplyr::between(log10fc, q25, q975)), 
    inherit.aes = F,
    aes(x = contrast, y = log10fc),
    alpha = 0.5,
    size = 0.5,
    width = 0.05
  ) +
  geom_errorbar(aes(ymin = q25, ymax = q975), width = 0.25, alpha = 0.75) +
  geom_point(shape = 21, alpha = 0.9, size = 3.5, aes(fill = signif)) +
  scale_fill_manual(values = c("grey", "#FEC368", "#FB8A30", "#F45C10")) +
  stat_pvalue_manual(
    data = sig %>% filter(sig != "ns"),
    label = "sig",
    xmin = "group1",
    xmax = "group2",
    y.position = "max",
    step.increase = 0.075,
    step.group.by = "Nuc",
    coord.flip = F,
    bracket.size = 0.3,
    tip.length = 0.02,
    label.size = 3,
    bracket.nudge.y = 0.75
  ) +
  theme_bw() + 
  facet_wrap(~Nuc, ncol = length(unique(nuc_comp$Nuc)), strip.position = "top") + 
  labs(
    title = "Nucleic acid type",
    x = "",
    y = "",
    size = "Num. subspecies",
    fill = "Sig. level \n[method/unconc]"  
  ) +
  theme(
    #plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
    legend.position = "none"
  )

# 4. Hypothesis testing for within method variation ----
plot_per_method <- rbind(
  plot_c2 %>% select(term, contrast, "type_val" = sizeclass, estimate, p.value, signif, conf.low, conf.high, q975, q25) %>% mutate(type = "size"),
  plot_e2 %>% select(term, contrast, "type_val" = envelope, estimate, p.value, signif, conf.low, conf.high,q975, q25) %>% mutate(type = "envelope"),
  plot_n2 %>% select(term, contrast, "type_val" = Nuc, estimate, p.value, signif, q975, conf.low, conf.high, q25) %>% mutate(type = "nucleic_acid"))

## Extract pairwise methods significance tests for each sizeclass level
sh <- avg_comparisons(lm2, variables = "Method", by = "sizeclass", hypothesis = "pairwise") %>% 
  separate(term, into = c("term1", "term2"), sep = "\\) - \\(")

sh2 <- sh %>% filter( # take only comparisons to the unconcentrated
  grepl("mean\\(MEM", term1) & grepl("mean\\(MEM", term2) | 
    grepl("mean\\(PEG", term1) & grepl("mean\\(PEG", term2) | 
    grepl("mean\\(VDC", term1) & grepl("mean\\(VDC", term2) | 
    grepl("mean\\(UF", term1) & grepl("mean\\(UF", term2) |   
    grepl("mean\\(NT", term1) & grepl("mean\\(NT", term2)
)

adj.p <- p.adjust(sh2$p.value, method = "holm")
sh2$p.value.adj <- adj.p
sh3 <- sh2 %>% filter(p.value.adj < 0.05) %>% select(term1, term2, estimate, p.value.adj)

## Envelope
eh <- avg_comparisons(lm2, variables = "Method", by = "envelope", hypothesis = "pairwise") %>% 
  separate(term, into = c("term1", "term2"), sep = "\\) - \\(")

eh2 <- eh %>% filter( # take only comparisons to the unconcentrated
  grepl("mean\\(MEM", term1) & grepl("mean\\(MEM", term2) | 
    grepl("mean\\(PEG", term1) & grepl("mean\\(PEG", term2) | 
    grepl("mean\\(VDC", term1) & grepl("mean\\(VDC", term2) | 
    grepl("mean\\(UF", term1) & grepl("mean\\(UF", term2) |   
    grepl("mean\\(NT", term1) & grepl("mean\\(NT", term2))

adj.p <- p.adjust(eh2$p.value, method = "holm")
eh2$p.value.adj <- adj.p
eh3 <- eh2 %>% filter(p.value.adj < 0.05) %>% select(term1, term2, estimate, p.value.adj)

## Nucleic acid
nh <- avg_comparisons(lm2, variables = "Method", by = "Nuc", hypothesis = "pairwise") %>% 
  separate(term, into = c("term1", "term2"), sep = "\\) - \\(")

nh2 <- nh %>% filter( # take only comparisons to the unconcentrated
  grepl("mean\\(MEM", term1) & grepl("mean\\(MEM", term2) | 
    grepl("mean\\(PEG", term1) & grepl("mean\\(PEG", term2) | 
    grepl("mean\\(VDC", term1) & grepl("mean\\(VDC", term2) | 
    grepl("mean\\(UF", term1) & grepl("mean\\(UF", term2) |   
    grepl("mean\\(NT", term1) & grepl("mean\\(NT", term2))

adj.p <- p.adjust(nh2$p.value, method = "holm")
nh2$p.value.adj <- adj.p
nh3 <- nh2 %>% filter(p.value.adj < 0.05) %>% select(term1, term2, estimate, p.value.adj)

# combine all adjusted pvalues
total_method_hyp <- rbind(sh3, eh3, nh3) %>% 
  separate(term1, into = c("term1", "group1"), sep = ", ") %>%
  separate(term2, into = c("term2", "group2"), sep = ", ") %>%
  mutate(group2 = str_remove(group2, "\\)")) %>%
  mutate(class_comparison = paste(group1, " - ", group2)) %>%
  mutate(contrast = str_remove_all(term1, "\\(|mean\\(|\\)")) %>%
  mutate(sig = case_when(
    p.value.adj <= 0.05 & p.value.adj >= 0.01 ~ "*",
    p.value.adj < 0.01 & p.value.adj >= 0.001 ~ "**",
    p.value.adj < 0.001 ~ "***",
    p.value.adj > 0.05 ~ "ns"))

max_val <- plot_per_method %>% ungroup() %>% summarize(max = max(estimate))# %>% dplyr::rename(group1 = type_val)
total_method_hyp <- total_method_hyp %>% mutate(max = max_val$max)

# Get jitter points from raw data
pdf_type <- pdf %>%
  select(Method, Replicate, subspecies, log10fc, sizeclass, envelope, Nuc) %>%
  pivot_longer(cols = c(sizeclass, envelope, Nuc), values_to = "type_val", names_to = "type") %>%
  mutate(contrast = paste0(Method, " - UNC")) %>%
  left_join(plot_per_method %>% select(contrast, type_val, q25, q975))

# Calculate coefficient of varition
coefv <- plot_per_method %>%
  group_by(contrast) %>%
  summarize(mean = mean(10^estimate), # transform back to original scale
            sd = sd(10^estimate),
            median = median(10^estimate),
            cv = sd / mean)

# Set plot factor order
plot_per_method$type_val <- factor(plot_per_method$type_val, levels = c("DNA", "RNA", "enveloped", "non-enveloped", "small", "mid", "large"))
plot_per_method$contrast <- factor(plot_per_method$contrast, levels = c("PEG - UNC", "VDC - UNC", "UF - UNC", "NT - UNC", "MEM - UNC"))

total_method_hyp$group1 <- factor(total_method_hyp$group1, levels = levels(plot_per_method$type_val))
total_method_hyp$group2 <- factor(total_method_hyp$group2, levels = levels(plot_per_method$type_val))

plot_per_method <- plot_per_method %>%
  mutate(signif = case_when(
    p.value > 0.05 ~ "ns",
    p.value <= 0.001 & estimate < 0 ~ "Sig-***",
    p.value <= 0.001 & estimate > 0 ~ "Sig+***",
    p.value <= 0.01 & estimate < 0 ~ "Sig-**",
    p.value <= 0.01 & estimate > 0 ~ "Sig+**",
    p.value <= 0.05 & estimate < 0 ~ "Sig-*",
    p.value <= 0.05 & estimate > 0 ~ "Sig+*"))

# Sign colors
sig_colors <- c(
  "ns" = "grey",
  "Sig+***" = "#F45C10",
  "Sig+**" = "#FB8A30",
  "Sig+*" = "#FEC368",
  "Sig-*" = "#A3C1E0",
  "Sig-*" = "#78AEDD",
  "Sig-***" = "#4D94DB"
)

### Plot with facet per method ----
plot_method <- ggplot(plot_per_method, aes(x = type_val, y = 10^estimate)) +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
  geom_hline(
    data = coefv, 
    aes(yintercept = median), 
    color = "skyblue3"
  ) +
  # geom_violin(
  #   data = pdf_type %>% group_by(contrast, type_val) %>% filter(dplyr::between(log10fc, q25, q975)), 
  #   aes(x = type_val, y = 10^log10fc),
  #   fill = "grey",
  #   color = "grey50",
  #   alpha = 0.7,
  # ) +
  geom_jitter(
    data = pdf_type, #%>% group_by(contrast, type_val) %>% filter(!dplyr::between(log10fc, q25, q975)), 
    inherit.aes = F,
    aes(x = type_val, y = 10^log10fc),
    alpha = 0.5,
    size = 0.5,
    width = 0.05
  ) +
  geom_errorbar(
    aes(ymin = 10^conf.low, ymax = 10^conf.high), 
    width = 0.25, 
    alpha = 0.75
  ) +
  geom_point(
    shape = 21, 
    alpha = 0.9, 
    size = 3.5, 
    aes(fill = signif)
  ) +
  stat_pvalue_manual(
    data = total_method_hyp,
    label = "sig",
    xmin = "group1",
    xmax = "group2",
    y.position = "max",
    step.increase = 0.05,
    step.group.by = "contrast",
    coord.flip = F,
    bracket.size = 0.3,
    tip.length = 0.02,
    label.size = 3,
    bracket.nudge.y = 1.1
  ) +
  geom_label(
    inherit.aes = F,
    data = coefv,
    aes(label = paste0(round(mean, 2), "±", round(sd, 2))),
    x = Inf,
    y = Inf,
    fill = "skyblue3",
    color = "white",
    vjust = 1.5,
    hjust = 1.1,
    size = 2.3
  ) +
  facet_wrap(~factor(contrast, levels = levels(plot_per_method$contrast)), ncol = 5) +
  labs(
    title = "Within method effects",
    x = "",
    y = "Estimate [fold change]",
    size = "Num. subspecies",
    fill = "Sig. level \n[method/unconc]"  ,
    #caption = "Blue line is median coefficient value. Dashed line is a fold change of 1 (eg. no enrichment)."
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = sig_colors) +
  scale_y_continuous(trans = "log10", labels = scales::label_number(), breaks = c(0.1, 1, 10, 100, 1000)) +
  scale_x_discrete(limits = levels(plot_per_method$type_val))

plot_method + canvas(height = 5, width = 12)