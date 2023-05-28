from selection_methods.gene_selection_shared import *
import time
import triku
import pandas as pd
# from scipy.sparse import issparse
import scanpy as sc


def preprocess_adata_triku(adata, size_factors="size_factors"):
    """Preprocess adata for method Triku

    Arguments
    ----------
    adata: AnnData
        adata object with raw counts in adata.X and normalisation size factors in adata.obs, reduced to highly variable genes
    size_factors: str
        Size factors key
    """

    sc.pp.filter_genes(adata, min_cells=5)

    # normalise(adata, key=size_factors)
    sc.pp.log1p(adata)

    return adata


def select_genes_triku(n, adata, size_factors="size_factors", **kwargs):
    """

    Parameters
    ----------
    n: int
        Number of genes to select, If None, the number is chosen automatically
    adata: AnnData
        adata object with normalized and logarithmized counts in adata.X
    size_factors: String (default: "size_factors")
        size_factors for normalization
    kwargs:
        use_raw : bool
            If True, selects the adata.raw, if it exists.
            To set the .raw propety, set as: adata.raw = adata.
            This matrix is adjusted to select the genes and cells that
            appear in the current adata. E.g. if we are running triku with a subpopulation, triku will select the cells
            from adata.raw of that subpopulation. If certain genes have been removed after saving the raw, triku will not
            consider the removed genes.
        n_divisions : int, None
            If the array of counts is not integer, number of bins in which each unit will be divided to account for
            that effect. For example, if n_divisions is 10, then 0.12 and 0.13 would be in the same bin, and 0.12 and 0.34
            in two different bins. If n_divisions is 2, the two cases would be in the same bin.
            The higher the number of divisions the more precise the calculation of distances will be. It will be more
            time consuming, though. If n_divisions is None, we will adjust it automatically.
        s : float
            Correction factor for automatic feature selection. Negative values imply a selction of more genes, and
            positive values imply a selection of fewer genes. We recommend values between -0.1 and 0.1.
        n_windows : int
            Number of windows used for median subtraction of Wasserstein distance.
        min_knn : int
            minimum number of expressed cells based on the knn to apply the convolution. If a gene has less than min_knn
            expressing cells, Wasserstein distance is set to 0, and the convolution is set as the knn expression.
        name: str
            Name of the run. If None, stores results in "triku_X". Else, stores it in "triku_X_{name}".
        verbose : str ['debug', 'triku', 'info', 'warning', 'error', 'critical']
            Logger verbosity output.

    """

    adata_cp = adata.copy()
    #adata_cp = preprocess_adata_triku(adata_cp, size_factors=size_factors)

    start = time.time()

    sc.pp.pca(adata_cp)
    sc.pp.neighbors(adata_cp, metric='cosine', n_neighbors=int(0.5 * len(adata_cp) ** 0.5))
    triku.tl.triku(adata_cp, n_features=n, **kwargs)

    selected_bool = adata_cp.var['highly_variable']

    end = time.time()
    took = end - start

    # get names of selected indices (Note: this was adjusted since we needed to add sc.pp.filter_genes in the preproc.)
    selectedGenes = adata_cp.var.index[selected_bool]
    selection = pd.DataFrame(index=adata.var_names,data={"idx": range(adata.shape[1]), "selection": adata.var_names})
    selection = selection.loc[selectedGenes]
    
    #adata.var['idx'] = range(adata.shape[1])
    #idx = adata.var['idx'][selected_bool]
    #selectedGenes = adata.var.index[selected_bool]
    #selection = pd.DataFrame({"idx": idx, "selection": selectedGenes})
    

    return selection, took
