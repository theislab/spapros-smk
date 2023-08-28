from selection_methods.gene_selection_shared import *
import scanpy as sc
import numpy as np
import pandas as pd
import time
import subprocess
import os

def preprocess_adata_scpnmf(adata, size_factors="size_factors", subsample=10000):
    """Preprocess adata for method scpnmf
    
    I actually just copied the preprocessing from scmer, we only included subsampling on top. Since scpnmf runs into 
    memory issues we use the same subsample size as for scmer.

    """

    normalise(adata, key=size_factors)
    sc.pp.log1p(adata)
    
    if subsample and (subsample < adata.n_obs):
        obs = np.random.choice(adata.n_obs, subsample, replace=False)
        adata = adata[obs]

    return adata

def select_genes_scpnmf(n, adata, r_exe, output_path, K=20, tol=1e-4, maxIter=1000, verboseN=True, zerotol=1e-10,
                        method="EucDist", label=None, mu=1, lam=0.01, seed=123,
                        toTest=True, cor_thres=0.7, pval_thres=0.01, return_fig=False, adj_method='BH', mc_cores=10,
                        ncol=4, toAnnotate=False, dim_use_basisSelect=None,
                        by_basis=False, return_trunW=True, dim_use_geneSelect=None, 
                        #conda_env="venv_scpnmf_R"
                        ):
    """Select genes via ScPNMF

        Arguments
        -----------
        n: int
            Number of genes to select
        adata: String
            adata object with raw counts in adata.X and path where an adata.h5ad is stored in uns["file_path"], NOTE: adata should not be normalized
        r_exe: String
            path to the Rscript.exe
        output_path:
            path where to store the output of the R script
        # for PNMF:
        K: int (default: 25)
            specification of the factorization rank (number of low dimentsions). They typically use a K between 10 and 20
        tol: float (default: 1e-3)
            a threshold bewlow which would be considered convergent
        maxIter: int (default: 500)
            number of max iteration times
        verboseN: bool (default: True)
            whether to print number of iterations
        zerotol: float (defualt: 1e-10)
            a threshold on basis loadings below which would be considered zero
        method: string (default: "EucDist")
            which method to be used. One of "EucDist", "KL" or "DPNMF"
        label: list (default: None)
            list whith cluster type for each cell, required only when method="DPNMF"
        mu: int (default: 1)
            controls the penalty term, larger mu represents heavier penalization of class distanced is DPNMF
        lam: fload (default: 0.01)
            contols the magnitude of within class distances, lager lam represents larger proportion of within calss distances in the total penaltry term
        seed: int (default: 123)
            random seed of the initialization
        # for basis selection:
        toTest: bool (default: True)
            whether to select bases by Person correlation w/ cell library size and test of multimodality
        cor_thres: float (default: 0.7)
            pearson correlation w/ cell library size cutoff
        pval_thres: float (default: 0.01)
            adjusted p-value cutoff on test of multimodality
        return_fig: bool (default: False)
            whether to print scatter plot of score vector against cell library size and distribution of score vectors for each basis
        adj_method: String (default: "BH")
            p-value correction method
        mc_cores: int (default: 10)
            the number of cores to use for function mclapply, for windows only 1 is supported!
        ncol: int (default: 4)
            columns for facets in plots
        toAnnotate: bool (default: False)
            whether to perform GO enrichment analysis on each basis
        dim_use_basisSelect: String (default: "NULL")
            comma-separated list of the bases (columns) to be used in the selected weight matrix, None value uses all bases
        # for getting informative genes
        by_basis: bool (default: False)
            return informative genes by basis or not
        return_trunW: bool (default: False)
            return the truncated weight matrix or not
        dim_use_geneSelect: String (default: "NULL")
            comma-separated list of the bases (columns) to be used in the selected weight matrix, None value uses all bases

    """
    # Get directory of this file
    file_path = os.path.realpath(__file__)
    dir_path = os.path.dirname(file_path)
    r_script = dir_path + "/gene_selection_scPNMF.R"

    start = time.time()
    command = r_exe + " " + r_script + " " + str(n) + " " + str(adata) + " " + str(K) + " " + str(
            tol) + " " + str(maxIter) + " " + str(verboseN) + " " + str(zerotol) + " " + str(method) + " " + str(
            label) + " " + str(mu) + " " + str(lam) + " " + str(seed) + " " + \
        str(toTest) + " " + str(cor_thres) + " " + str(pval_thres) + " " + str(return_fig) + " " + str(
            adj_method) + " " + str(mc_cores) + " " + str(ncol) + " " + str(toAnnotate) + " " + str(
            dim_use_basisSelect) + " " + \
        str(by_basis) + " " + str(return_trunW) + " " + str(dim_use_geneSelect) + " " + \
        output_path #+ " " + conda_env
    print(command)
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    process.wait()
    end = time.time()
    took = end - start
    adata = sc.read(adata)
    selectedGenes = list(pd.read_csv(output_path, sep="\t")["x"])
    adata.var['idx'] = range(adata.shape[1])
    idx = [int(adata.var['idx'][adata.var.index == x]) for x in selectedGenes]
    selection = pd.DataFrame({"idx": idx, "selection": selectedGenes})

    # example command
    # C:\Users\stras\anaconda3\envs\conda_env_gene_selection\lib\Rscript gene_selection_scPNMF.R 25 C:\Users\stras\Documents\8.Semester\Louis_Hiwi\gene_selection_methods\examples\data\pbmc3k_highly_variable4000.h5ad 20 0.0001 1000 True 1e-10 EucDist leiden 1 0.01 123 True 0.7 0.01 False BH 10 4 False None False True None C:\Users\stras\Documents\8.Semester\Louis_Hiwi\gene_selection_methods\examples\examples/selected_genes\pbmc3k_scpnmf_25_temp.tsv
    # /big/st/strasserl/gene_selection_methods/parent_dir/environments/venv_scpnmf_R/bin/Rscript selection_methods/gene_selection_scPNMF.R 10 None 20 0.0001 1000 True 1e-10 EucDist leiden 1 0.01 123 True 0.7 0.01 False BH 10 4 False None False True None selections/tmp/scpnmf_tmp.tsv

    return selection, took
