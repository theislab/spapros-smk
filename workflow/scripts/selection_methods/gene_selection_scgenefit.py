from scGeneFit.functions import *
from selection_methods.gene_selection_shared import *
import scanpy as sc
import numpy as np
import pandas as pd
import time


def preprocess_adata_scgenefit(adata, size_factors="size_factors"):
    """Preprocess adata for method ScGeneFit

    Arguments
    ---------
    adata: AnnData
        adata object with raw counts in adata.X and normalisation size factors in adata.obs, reduced to highly variable genes
    label: str
        Celltype clustering key
    hyper_parameter1: float
        Filters out genes that are not expressed in `hyper_parameter1` % cells of at least one
        cluster in groupbs of adata.obs[groupb_by]

    """

    normalise(adata, key=size_factors)
    sc.pp.log1p(adata)
    # xxxpy.filter_genes(adata,hyper_parameter1=0.75,group_by=group_by)

    # convert adata to Nxd numpy array where N = number of points (genes) and d = dimensions (cells)
    # data = adata.X.toarray()
    return adata


def select_genes_scgenefit(n, adata, label, hierarchical=False, size_factors="size_factors", method='centers', epsilon=1, sampling_rate=1,
                           n_neighbors=3, max_constraints=1000, redundancy=0.01, verbose=False):
    """Select genes via method ScGeneFit

    Arguments
    ---------
    n: int
        Number of genes to select
    adata: AnnData
        adata object with normalized and logarithmized counts in adata.X
    label: String or list
        name of adata.obs where a list of labels can be found: eg. cell type, cell state or
        list with labels (N labels, one per point (gene names))
        if hierarchy = True:  list with T lists of labels, where T is the number of layers in the hierarchy (N labels per list, one per point)
    hierarchical: bool
        whether the hierarchical method should be used. If True, the labels have to have a hierarchy
    hyper_parameter1:
        TODO
    method: String (default 'centers')
        'centers', 'pairwise', or 'pairwise_centers'
    epsilon: int (default: 1)
        constraints will be of the form expr>Delta, where Delta is chosen to be epsilon times the norm of the smallest constraint (default 1)
        (This is the most important parameter in this problem, it determines the scale of the constraints,
        the rest the rest of the parameters only determine the size of the LP)
    sampling_rate: float (default: 1)
        (if method=='pairwise' or 'pairwise_centers') selects constraints from a random sample of proportion sampling_rate (default 1)
    n_neighbors: int (default: 3)
        (if method=='pairwise') chooses the constraints from n_neighbors nearest neighbors (default 3)
    max_constraints: int (default: 1000)
        maximum number of constraints to consider (default 1000)
    redundancy: float (default: 0.01)
        (if method=='centers') in this case not all pairwise constraints are considered
        but just between centers of consecutive labels plus a random fraction of constraints given by redundancy
        if redundancy==1 all constraints between pairs of centers are considered
    verbose: bool (default: False)

    """

    #adata = preprocess_adata_scgenefit(adata, hyper_parameter1, size_factors=size_factors)

    np.random.seed(0) # As recommended in the scgenefit repo, however there seems to be additional randomness
                      # when running wih hierarchical = False.

    if type(label) == str:
        labels = np.array(adata.obs[label])

    # np.ndarray needed
    if not type(adata.X) == np.ndarray:
        adata_cp = adata.copy().X.toarray()
    else:
        adata_cp = adata.copy().X

    # hierarchical method
    if hierarchical:
        start = time.time()


        idx = get_markers_hierarchy(data=adata_cp, labels=labels, num_markers=n, method=method,
                                              epsilon=epsilon, sampling_rate=sampling_rate, n_neighbors=n_neighbors,
                                              max_constraints=max_constraints, redundancy=redundancy, verbose=verbose)
        end = time.time()
    # non-hierarchical variant
    else:
        start = time.time()
        idx = get_markers(data=adata_cp, labels=labels, num_markers=n, method=method,
                                    epsilon=epsilon, sampling_rate=sampling_rate, n_neighbors=n_neighbors,
                                    max_constraints=max_constraints, redundancy=redundancy, verbose=verbose)
        end = time.time()

    took = end - start

    # selected = pd.Series([False] * len(adata))
    # selected.iloc[selectedGenes] = True
    #
    # selection = pd.DataFrame(index=adata.var.index)
    # selection["selection"] = selected

    adata.var['idx'] = range(adata.shape[1])
    selectedGenes = adata.var.index[idx].values
    selection = pd.DataFrame({"idx": idx, "selection": selectedGenes})

    return selection, took



