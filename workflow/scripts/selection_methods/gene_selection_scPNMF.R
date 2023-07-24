# this is a script that can be invoked from the python wrapper to do gene selection with the R package scPNMF
# for infos on the arguments ets see docstrings of the python wrapper

# load libraries
suppressPackageStartupMessages(library(scPNMF))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(umap))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reticulate))


# parse arguments
args <- commandArgs(trailingOnly = TRUE)
n <- as.integer(args[1])
input_path <- args[2]
K  <- as.integer(args[3])
tol <- as.double(args[4])
maxIter <- as.double(args[5])
verboseN <- as.logical(args[6])
zerotol <- as.double(args[7])
method <- args[8]
label <- args[9]
mu <- as.double(args[10])
lam <- as.double(args[11])
seed <- as.integer(args[12])
toTest <- as.logical(args[13])
cor_thres <- as.double(args[14])
pval_thres <- as.double(args[15])
return_fig <- as.logical(args[16])
adj_method <- args[17]
mc_cores <- as.integer(args[18])
ncol <- as.integer(args[19])
toAnnotate <- as.logical(args[20])
dim_use_basisSelect <- args[21]
if(dim_use_basisSelect == "None"){
  dim_use_basisSelect <- NULL
}
by_basis <- as.logical(args[22])
return_trunW <- as.logical(args[23])
dim_use_geneSelect <- args[24]
if(dim_use_geneSelect == "None"){
  dim_use_geneSelect <- NULL
}
output_path <- args[25]
#conda_env <- args[26]

# specify conda env
#use_condaenv(condaenv = conda_env, required = T)
suppressPackageStartupMessages(library(anndata))

# read input (raw counts)
input <- read_h5ad(input_path)

# transpose matrix (cells x genes --> genes x cells), convert sparse to dense matrix and logarithmize the counts
log_matrix <- log1p(as.matrix(input$T$X))

# extract labels from anndata
labls <- input$obs[[label]]

# run PNMF
res_pnmf <- scPNMF::PNMFfun(X = log_matrix, K=K, tol=tol, maxIter=maxIter, verboseN=verboseN, zerotol=zerotol, method=method, label=labls, mu=mu, lambda=lam, seed=seed)

# extract weiths and scores from PNMF
W <- res_pnmf$Weight
S <- res_pnmf$Score

# select bases
print(dim_use_basisSelect)
W_select <- scPNMF::basisSelect(W = W, S = S, X=log_matrix, toTest=toTest, cor_thres=cor_thres, pval_thres=pval_thres, return_fig=return_fig, adj_method=adj_method, ncol=ncol, toAnnotate=toAnnotate, dim_use=dim_use_basisSelect)

# get the n most informative genes according to the selected bases
print(dim_use_geneSelect)
ig <- getInfoGene(W_select, M = n, by_basis=by_basis, return_trunW=return_trunW, dim_use=dim_use_geneSelect)

# write (tmp) output file
write.table(ig$InfoGene, file=output_path, quote=F, sep="\t")

# print(input_path)
# print(output_path)
# print(label)
# print(n)
# print(conda_env)


# minimal example
# n <- 25
# label <- "leiden"
# # input_path = "C:\\Users\\stras\\Documents\\8.Semester\\Louis_Hiwi\\gene_selection_methods\\examples\\data\\pbmc3k_highly_variable4000.h5ad"
# input <- read_h5ad(input_path)
# log_matrix <- log1p(as.matrix(input$T$X))
# labls <- input$obs[[label]]
# res_pnmf <- scPNMF::PNMFfun(X = log_matrix, K=10, maxIter=500, method="DPNMF", label=labls)
# W <- res_pnmf$Weight
# S <- res_pnmf$Score
# W_select <- scPNMF::basisSelect(W = W, S = S, X=log_matrix, mc.cores=1)
# colnames(W_select)
# ig <- getInfoGene(W_select, M = n)
# # output_path <- "C:\\Users\\stras\\Documents\\8.Semester\\Louis_Hiwi\\gene_selection_methods\\pbmc3k_scpnmf_25_DPNMF_temp.tsv"
# write.table(ig$InfoGene, file=output_path, quote=F, sep="\t")
