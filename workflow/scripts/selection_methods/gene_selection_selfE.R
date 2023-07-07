# this is a script that can be invoked from the python wrapper to do gene selection with the R package scPNMF
# for infos on the arguments ets see docstrings of the python wrapper

## Set libPaths for jobs on icb cluster
#if (.libPaths() == "/opt/R/lib/R/library") {
#    .libPaths(c("/home/louis.kuemmerle/bin",.libPaths()))
#}

# load libraries
suppressPackageStartupMessages(library(reticulate))
#if (.libPaths()[1] == "/home/louis.kuemmerle/bin") {
#    #use_python("/opt/python/bin/python")
#    #use_python("/opt/python/bin/python3")
#    #use_python("/opt/python/lib/python3.8")
#    use_python("/home/louis.kuemmerle/.local/share/r-miniconda/envs/r-reticulate/bin/python")
#}


# parse arguments
args <- commandArgs(trailingOnly = TRUE)
n <- as.integer(args[1])
prep_adata_path <-  args[2]
output_path <- args[3]
dir_of_script <- args[4]
#conda_env <- args[4]

## example:
# n <- 25
# prep_adata_path <- "selections/tmp/preprocessed_adata.h5ad"
# output_path <- "selections/tmp/selfe_tmp.tsv"
# conda_env <- "/big/st/strasserl/gene_selection_methods/gene_selection_methods/environments/venv_selfe"


## specify conda env
#use_condaenv(condaenv = conda_env, required = T)
#
# load or install packages that are not available from conda thus not installed yet
packages <- c("anndata", "pracma")
lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

# read input (raw counts), samples are in row and features are in column
input <- read_h5ad(prep_adata_path)

# load gene selection method
source(paste0(dir_of_script,"/SelfE.R"))

# selected_genes: indices one-based
selected_genes <- SelfE(input$X, n)

# selected_genes_0: indices zero-based
selected_genes_0 <- selected_genes-1

write.table(selected_genes_0, file=output_path, quote=F, sep="\t", row.names=F)

