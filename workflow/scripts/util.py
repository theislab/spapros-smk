import warnings
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse, csr_matrix
from sklearn.utils import sparsefuncs
from typing import Optional, List

def normalise(adata, key='size_factors'):
    """Normalise raw counts adata with size factors

    adata: AnnData
        Contains size factors in adata.obs[key]
    key: str
        Key for size factors
    """

    if issparse(adata.X):
        sparsefuncs.inplace_row_scale(adata.X, 1 / adata.obs[key].values)
    else:
        adata.X /= adata.obs[key].values[:, None]


def preprocess_adata(adata, size_factors="size_factors"):
    """Preprocess adata for method SCMER

    Arguments
    ----------
    adata: AnnData
        adata object with raw counts in adata.X and normalisation size factors in adata.obs, reduced to highly variable genes
    size_factors: str
        Size factors key
    """


    normalise(adata, key=size_factors)
    sc.pp.log1p(adata)

    return adata


def get_highly_variable(adata,N=8000,key_added='highly_variable'):
    """Compute N highly variable genes of given dataset for given normalisation
    
    A copy of the raw data is normalised and log1p transformed to extract N highly variable genes.
    """
    a = adata.copy()
    a.X /= a.obs['size_factors'].values[:,None]
    a.X = csr_matrix(a.X)
    sc.pp.log1p(a)
    sc.pp.highly_variable_genes(a, n_top_genes=N, flavor='cell_ranger', inplace=True)
    adata.var[key_added] = a.var['highly_variable']


def subset_to_n_celltypes(adata, n : int, ct_key : str = "celltype", exact : bool = False):
    """Subset adata to maximally n cell types, cell types with more cells are chosen first
    
    Arguments
    ----------
    adata: AnnData
        adata object
    n: int
        Number of cell types to select
    ct_key: str
        Key for cell type annotation
    exact: bool
        Whether to raise an error if there are less than n cell types in adata
        
    Returns
    -------
    adata: AnnData
        adata object with only n cell types
    """
    if ct_key is None:
        raise ValueError("ct_key must be specified when subsetting adata to n cell types")
    
    ct_counts = adata.obs[ct_key].value_counts()
    if (len(ct_counts) < n) and exact:
        raise ValueError(f"adata contains only {len(ct_counts)} cell types, but {n} were requested")
    elif len(ct_counts) < n:
        print(f"adata contains only {len(ct_counts)} cell types, but {n} were requested. Moving on with all cell types")
        n = len(ct_counts)
    ct_counts = ct_counts.iloc[:n]
    ct_names = ct_counts.index.values
    adata = adata[adata.obs[ct_key].isin(ct_names), :]
    
    return adata


def sample_n_cells_per_ct(adata, n, ct_key="celltype", seed=0):
    """Sample n cells per cell type (keep all cells if there are less than n)
    
    Arguments
    ---------
    adata: AnnData
        adata object
    n: int
        Number of cells to sample per cell type
    ct_key: str
        Key for cell type annotation
    seed: int
        Random seed
        
    Returns
    -------
    adata: AnnData
        adata object with n cells per cell type (or all cells if there are less than n)
    """
    
    if ct_key is None:
        raise ValueError("ct_key must be specified when sampling n cells per cell type")
    
    ct_counts = adata.obs[ct_key].value_counts()
    cts = ct_counts.index.values
    
    obs = []
    for ct in cts:
        if ct_counts[ct] > n:
            np.random.seed(seed)
            ct_obs = np.random.choice(adata.obs_names[adata.obs[ct_key] == ct], n, replace=False)
            obs += list(ct_obs)
        else:
            obs += list(adata.obs_names[adata.obs[ct_key] == ct])
    
    adata = adata[obs, :]
    
    return adata

def get_selection_df(adata, selected_genes: list) -> pd.DataFrame:
    """Get selection dataframe
    
    selection scripts for external methods have this output, let's stick to it.
    
    """
    adata.var['idx'] = range(adata.shape[1])
    idx = [int(adata.var['idx'][adata.var.index == x]) for x in selected_genes]
    selection = pd.DataFrame({"idx": idx, "selection": selected_genes})
    
    return selection


def bootstrap_sample(adata, obs_key=None, noise_level=1e-3, seed=None, obs_names_unique=True):
    """
    Generate a bootstrap sample from an AnnData object. If obs_key is provided, 
    preserves the category frequencies of obs_key. Introduces small noise to duplicated samples.
    
    Arguments
    ---------
    Parameters:
    adata: AnnData
        The annotated data matrix.
    obs_key: str, optional
        Key for the observations to keep category frequencies constant. If None, the whole dataset is bootstrapped 
        without preserving categories.
    noise_level: float, optional
        The magnitude of noise to introduce to duplicated samples. Defaults to 1e-3.
    seed: int, optional
        Random seed for reproducibility. Defaults to None.
    obs_names_unique: bool, optional
        Whether to make observation names unique. Defaults to True.
    
    Returns
    -------
    AnnData: A bootstrap sample of the AnnData object.
    """
    np.random.seed(seed)
    
    if obs_key:
        sampled_indices = []

        # Loop through each unique category in obs_key
        for category in adata.obs[obs_key].unique():
            category_indices = adata.obs.index[adata.obs[obs_key] == category].tolist()
            bootstrap_indices = np.random.choice(category_indices, size=len(category_indices), replace=True)
            sampled_indices.extend(bootstrap_indices)
    else:
        sampled_indices = np.random.choice(adata.obs.index, size=adata.n_obs, replace=True)
    
    with warnings.catch_warnings():
        if obs_names_unique:
            warnings.filterwarnings("ignore", category=UserWarning, module="anndata._core.anndata")
        bootstrap_adata = adata[sampled_indices, :].copy()
    
    # Identify duplicates
    _, idx, counts = np.unique(bootstrap_adata.obs_names, return_index=True, return_counts=True)
    duplicate_indices = idx[counts > 1]
    
    # Convert to dense if necessary
    is_sparse = issparse(bootstrap_adata.X)
    if is_sparse:
        bootstrap_adata.X = bootstrap_adata.X.toarray()
    
    # Add noise
    noise = np.random.normal(0, noise_level, bootstrap_adata.X.shape)
    noise[bootstrap_adata.X == 0] = 0
    bootstrap_adata.X[duplicate_indices, :] += noise[duplicate_indices, :]
    
    # Convert back to sparse if necessary
    if is_sparse:
        bootstrap_adata.X = csr_matrix(bootstrap_adata.X)
    
    # Make observation names unique
    if obs_names_unique:
        bootstrap_adata.obs_names_make_unique(join="-")
    
    return bootstrap_adata


def clean_adata(
    adata: sc.AnnData,
    obs_keys: Optional[List[str]] = None,
    var_keys: Optional[List[str]] = None,
    uns_keys: Optional[List[str]] = None,
    obsm_keys: Optional[List[str]] = None,
    varm_keys: Optional[List[str]] = None,
    obsp_keys: Optional[List[str]] = None,
    inplace: bool = True,
) -> Optional[sc.AnnData]:
    """Removes unwanted attributes of an adata object.

    Args:
        adata:
            Anndata object.
        obs_keys:
            Columns in `adata.obs` to keep.
        var_keys:
            Columns in `adata.var` to keep.
        uns_keys:
            Keys of `adata.uns` to keep.
        obsm_keys:
            Columns of `adata.obsm` to keep.
        varm_keys:
            Columns of `adata.varm` to keep.
        obsp_keys:
            Keys of `adata.obsp` to keep.
        inplace:

    Returns:
        If not inpace, the cleaned Anndata without any annotation.
    """
    if not obs_keys:
        obs_keys = []
    if not var_keys:
        var_keys = []
    if not uns_keys:
        uns_keys = []
    if not obsm_keys:
        obsm_keys = []
    if not varm_keys:
        varm_keys = []

    if inplace:
        a = adata
    else:
        a = adata.copy()
    for obs in [o for o in a.obs_keys() if o not in obs_keys]:
        del a.obs[obs]
    for var in [v for v in a.var_keys() if v not in var_keys]:
        del a.var[var]
    for uns in [u for u in a.uns_keys() if u not in uns_keys]:
        del a.uns[uns]
    for obsm in [om for om in a.obsm_keys() if om not in obsm_keys]:
        del a.obsm[obsm]
    for varm in [vm for vm in a.varm_keys() if vm not in varm_keys]:
        del a.varm[varm]
    # for obsp in [op for op in a.obsp_keys() if not op in obsp_keys]:
    # 	del a.obsp[obsp]
    if not inplace:
        return a
    else:
        return None



def select_pca_genes(
    adata: sc.AnnData,
    n: int,
    absolute: bool = True,
    n_pcs: int = 20,
    inplace: bool = True,
) -> pd.DataFrame:
    """Select n features based on pca loadings.

    Args:
        adata:
            Data with log normalised counts in adata.X.
        n:
            Number of selected features.
        absolute:
            Take absolute value of loadings.
        n_pcs:
            Number of PCs used to calculate loadings sums.
        inplace:
            Save results in :attr:`adata.var` or return dataframe.

    Returns:

        - if not inplace:
            pd.DataFrame (like adata.var) with columns:

                - 'selection': bool indicator of selected genes
                - 'selection_score': pca loadings based score of each gene
                - 'selection_ranking': ranking according selection scores

        - if inplace:
            Save results in `adata.var[['selection','selection_score','selection_ranking']]`.

    """

    a = adata.copy()

    if n_pcs > a.n_vars:
        n_pcs = a.n_vars

    clean_adata(a)

    sc.pp.pca(
        a,
        n_comps=n_pcs,
        zero_center=True,
        svd_solver="arpack",
        random_state=0,
        return_info=True,
        copy=False,
    )

    loadings = a.varm["PCs"].copy()[:, :n_pcs]
    if absolute:
        loadings = abs(loadings)

    scores = pd.DataFrame(index=adata.var.index, data={"scores": np.sum(loadings, axis=1)})

    selected_genes = scores.nlargest(n, "scores").index.values
    selection = pd.DataFrame(
        index=scores.index,
        data={
            "selection": False,
            "selection_score": scores["scores"],
            "selection_ranking": scores["scores"].rank(method="dense", ascending=False),
        },
    )
    selection.loc[selected_genes, "selection"] = True

    if inplace:
        adata.var[["selection", "selection_score", "selection_ranking"]] = selection[
            ["selection", "selection_score", "selection_ranking"]
        ]
    else:
        return selection
