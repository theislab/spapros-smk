import scanpy as sc
from scipy.sparse import issparse
from sklearn.utils import sparsefuncs


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