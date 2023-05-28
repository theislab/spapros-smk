# wrapper for gene selection with geneBasis

# load libraries
#suppressPackageStartupMessages(library(devtools))
#install_github('MarioniLab/geneBasisR')
suppressPackageStartupMessages(library(geneBasisR))
suppressPackageStartupMessages(library(zellkonverter))
suppressPackageStartupMessages(library(SingleCellExperiment))

# parse arguments

args <- commandArgs(trailingOnly = TRUE)
n <- as.integer(args[1])
prep_adata_path <-  args[2]
output_path <- args[3]

if(args[4] == 'None'){
  genes_base <- NULL
}else{
  genes_base <- strsplit(args[4], split = " ")
}
batch <- args[5]
if(batch == "None"){
  batch <- NULL
}
neigh <- as.integer(args[6])
if(is.na(neigh)){
  neigh <- 5
}
mink <- as.integer(args[7])
if(mink == "None"){
  mink <- 3
}
npc_selection <- as.integer(args[8])
if(is.na(npc_selection)){
  npc_selection <- NULL
}
npc_all <- as.integer(args[9])
if(is.na(npc_all)){
  npc_all <- 50
}
if(args[10] == 'None'){
  genes_discard <- NULL
}else{
  genes_discard <- strsplit(args[10], split = " ")
}
if(args[11] == 'None'){
  genes_discard_prefix <- NULL
}else{
  genes_discard_prefix <- strsplit(args[11], split = " ")
}
verbose <- as.logical(args[12])
if(verbose == "None"){
  verbose <- TRUE
}
#conda_env <- args[13]

stat_all <- NULL
#' @param stat_all If True graph and corresponding Minkowski distances have been calculated prior to search, provide this data here.
#' It can be useful if gene_search is desired to be recycled (e.g. for selecting multiple libraries with different inputs such as n_genes_total and genes_base)
#' Ensure that colnames = c("gene", "dist_all"). Default stat_all=NULL - in case this info is not supplied.


# minimal example:
# conda_env <- "/big/st/strasserl/gene_selection_methods/parent_dir/environments/venv_genebasis"
# prep_adata_path <- "/big/st/strasserl/gene_selection_methods/parent_dir/selections/tmp/preprocessed_adata.h5ad"
# output_path <- "/big/st/strasserl/gene_selection_methods/parent_dir/selections/tmp/genebasis_tmp.tsv"

## specify conda env
suppressPackageStartupMessages(library(reticulate))
#use_condaenv(condaenv = conda_env, required = T)
use_python("/usr/bin/python")

# read input (log normalized HVGs)
# creates a new conda env in ~/.cache/basilisk/1.2.1/zellkonverter-1.0.3/anndata_env
input <- readH5AD(prep_adata_path)   # basilisk
# rows: genes, columns: cells, $cell: unqique cell ids, assay called logcounts
logcounts(input) <- assay(input)
input$cell <- colnames(input)

#' n_row = 1000
#' n_col = 100
#' sce = SingleCellExperiment(assays = list(logcounts = matrix(rnorm(n_row*n_col), ncol=n_col)))
#' rownames(sce) = as.factor(1:n_row)
#' colnames(sce) = c(1:n_col)
#' sce$cell = colnames(sce)
#' genes = rownames(sce)
#' out = gene_search(sce, n_genes_total = 5)
#'
#' SingleCellExperiment object, wehre normalized counts stored in 'logcounts' assay

selection <- gene_search(input,
                        n_genes_total = n,
                        genes_base = genes_base,
                        batch = batch,
                        n.neigh = neigh,
                        p.minkowski = mink,
                        nPC.selection = npc_selection,
                        nPC.all = npc_all,
                        genes.discard = genes_discard,
                        genes.discard_prefix = genes_discard_prefix,
                        verbose = verbose,
                        stat_all = stat_all)

print(selection)

# write (tmp) output file
write.table(selection, file=output_path, quote=F, sep="\t", row.names=F)