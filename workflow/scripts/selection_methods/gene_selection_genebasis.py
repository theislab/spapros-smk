import pandas as pd
import numpy as np
import scanpy as sc
import time
import subprocess
import os

from selection_methods.gene_selection_shared import *


def select_genes_genebasis(n, adata, tmp_dir, #conda_env, 
                           r_exe="Rscript", #size_factors="size_factors",
                           genes_base=None, batch=None, neigh=5, mink=3, npc_selection=None, npc_all=50,
                           genes_discard=None, genes_discard_prefix=None, verbose=True):
    """

    Parameters
    ----------
    n: int
        Number of genes to select
    adata: Anndata
        Anndata object with raw counts in adata.X and size factors in adata.obs[size_factors]
    tmp_dir: String
        path where tmp files should be stored
    conda_env:
        name (or prefix) of conda environment
    r_exe:
        path to Rscript (in conda env)
    size_factors: String (default: "size_factors")
        name of adata.obs were size factors are stored
    genes_base: list of Strings (default: None)
        Character vector specifying base genes to construct first Selection graph. Default=NULL in case no genes are supplied.
    batch: String (default: None)
        Name of the field in colData(sce) to specify batch.
    neigh: int (default: 5)
        Positive integer > 1, specifying number of neighbors to use for kNN-graph.
    mink: int (default: 3)
        p.minkowski Order of Minkowski distance.
    npc_selection: int (default: None)
        Scalar specifying number of PCs to use for Selection Graphs.
    npc_all: int (default: 50)
        Scalar specifying number of PCs to use for True Graph.
    genes_discard: list of Strings (default: None)
        Character vector containing genes to be excluded from candidates (note that they still will be used for graphs construction. If you want to exclude them from graph construction as well, just discard them prior in sce object). Default = NULL and no genes will be discarded.
    genes_discard_prefix: list of Strings (default: None)
        Character vector containing prefixes of genes to be excluded (e.g. Rpl for L ribosomal proteins. Note that they still will be used for graphs construction. If you want to exclude them from graph construction as well, just discard them prior in sce object). Default = NULL and no genes will be discarded.
    verbose: Boolean (default: True)
        Boolean identifying whether intermediate print outputs should be returned
    """

    #adata = preprocess_adata(adata, size_factors=size_factors)
    preprocessed_adata_path = os.path.join(tmp_dir, "preprocessed_adata.h5ad")
    adata.write_h5ad(preprocessed_adata_path)
    del adata  # need to delete the anndata because otherwise the h5ad file is locked and can't be read by the R script

    select_scripts_dir = os.path.dirname(os.path.abspath(__file__))

    start = time.time()
    command = r_exe + f" {select_scripts_dir}/gene_selection_geneBasis.R " + " ".join([str(n),  # 1
                                                                                   preprocessed_adata_path,  # 2
                                                                                   os.path.join(tmp_dir, "genebasis_tmp.tsv"),  # 3
                                                                                   str(" ".join(genes_base)) if genes_base is not None else str(None),  # 4
                                                                                   str(batch),  # 5
                                                                                   str(neigh),  # 6
                                                                                   str(mink),  # 7
                                                                                   str(npc_selection),  # 8
                                                                                   str(npc_all),  # 9
                                                                                   str(" ".join(genes_discard)) if genes_discard is not None else str(None),  # 10
                                                                                   str(" ".join(genes_discard_prefix)) if genes_discard_prefix is not None else str(None),  # 11
                                                                                   str(verbose)#,  # 12
                                                                                   #str(conda_env)  # 13
                                                                                   ])

    # example command:
    # ./environments/venv_genebasis/bin/Rscript selection_methods/gene_selection_geneBasis.R 25 selections/tmp/preprocessed_adata.h5ad selections/tmp None None 5 3 None 50 None None True

    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    process.wait()
    end = time.time()
    took = end - start
    adata = sc.read(preprocessed_adata_path)
    selectedGenes = list(pd.read_csv(os.path.join(tmp_dir,"genebasis_tmp.tsv"), sep="\t")['gene'])
    adata.var['idx'] = range(adata.shape[1])
    idx = [int(adata.var['idx'][adata.var.index == x]) for x in selectedGenes]
    selection = pd.DataFrame({"idx": idx, "selection": selectedGenes})


    return selection, took