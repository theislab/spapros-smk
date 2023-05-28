import scanpy as sc
from scipy.sparse import issparse
from sklearn.utils import sparsefuncs
import numpy as np
from sklearn.neighbors import NearestCentroid


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


def performance(n, adata, markers, label="leiden"):
    if issparse(adata.X):
        data = adata.copy().X.toarray()
    else:
        data = adata.copy().X
    clf = NearestCentroid()
    labels = np.array(adata.obs[label])

    clf.fit(data, labels)
    accuracy = clf.score(data, labels)
    print("Accuracy (whole data,", data.shape[1], " markers): ", accuracy)

    clf.fit(data[:, markers], labels)
    accuracy_markers = clf.score(data[:, markers], labels)
    print("Accuracy (selected", n, "markers)", accuracy_markers)






