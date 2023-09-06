# this script is for building source tree for R packages
# author: laojp
# time: 2023.09.06
# position: SYSUCC bioinformatic platform
# uasge 
#  Rscript build_source_tree.R [seed_file] [output_dir]

setwd(commandArgs()[2])
print("-------- firct check the packages exist or not --------")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  print("BiocManager is not installed and install")
  install.packages("BiocManager", repos = NULL)
}
for (i in c("Biostrings","BSgenome")) {
  if (!require(i, character.only = TRUE)) {
    BiocManager::install(i, ask = FALSE)
  }
}

print("-------- make source tree for packages --------")
forgeBSgenomeDataPkg(commandArgs()[1])
