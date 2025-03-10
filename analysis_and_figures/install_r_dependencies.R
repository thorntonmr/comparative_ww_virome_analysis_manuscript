# Install all libraries needed to produce figure draft
required_packages <- c(
  "conflicted", 
  "tidyverse", 
  "Biostrings", 
  "reactable", 
  "reactablefmtr", 
  "viridis", 
  "data.table"
  "ggview",
  "ggthemes",
  "ggrepel",
  "ggpubr",
  "viridis",
  "patchwork",
  "MASS",
  "performance",
  "marginaleffects",
  "data.table",
  "vegan",
  "rio",
  "RVAideMemoire",
  "lmtest",
  "knitr"
  )

cat("Checking package installations...\n")
for (package_name in required_packages) {
  if (package_name %in% installed.packages()) {
    cat("  - Package", package_name, "is installed. \n")
  } else{
    cat("WARNING: The package", package_name, "is not installed. Installing now...\n")
    suppressMessages(install.packages(package_name, dependencies = TRUE))   
  }
}

cat("Checking packages for report creation...\n")
if(!require(dataui)){
  cat("R: WARNING: Could not find dataui package. Installing now...\n")
  remotes::install_github("timelyportfolio/dataui")
}
if(!require(bedtoolsr)){
  cat("R: WARNING: Could not find rbedtools package. Installing now...\n")
  devtools::install_github("PhanstielLab/bedtoolsr")
}

cat("All packages successfully installed.")
