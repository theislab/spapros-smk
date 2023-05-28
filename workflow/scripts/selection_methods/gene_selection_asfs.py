from selection_methods.ASFS_min_complexity import SVM_active_feature_selection as SVM_active_feature_selection_min_complexity
from selection_methods.ASFS_min_cell import SVM_active_feature_selection as SVM_active_feature_selection_min_cell
from selection_methods.gene_selection_shared import *
import time
import numpy as np
import pandas as pd
import random
import sys


def select_genes_asfs(n, adata, label, num_samples, balance=None, method="min_cell", size_factors="size_factors"):
    """

    Parameters
    ----------
    n: int
        Number of genes to select
    adata: AnnData
        adata object with normalized and logarithmized counts in adata.X
    label: String
        name of column in adata.obs, where cell types clusters are stored
    balance: bool
        balance the number of cells of each class or just randomly select cells at each loop (only needed when method=="min_complexity")
    num_samples: int
        the number of cells we would use at each loop
    method: String
        either "min_complexity" or "min_cell"
    size_factors: String
        size_factors for normalization
    """

    #adata = preprocess_adata(adata, size_factors=size_factors)
    if not type(adata.X) == np.ndarray:
        data = adata.copy().X.toarray()
    else:
        data = adata.copy().X

    #target = adata.obs[label].values.reshape((data.shape[0],)).astype(np.uint8) - 1
    target = adata.obs[label].astype('category').cat.codes.values.astype(np.uint8) - 1    

    idx = np.arange(np.shape(data)[0])
    random.shuffle(idx)
    X_train = data[idx[:int(np.shape(data)[0] * 4 / 5)], :]
    y_train = target[idx[:int(np.shape(data)[0] * 4 / 5)]]
    X_test = data[idx[int(np.shape(data)[0] * 4 / 5):], :]
    y_test = target[idx[int(np.shape(data)[0] * 4 / 5):]]

    start = time.time()

    if method == "min_complexity":
        idx, num_samples_list, train_errors, test_errors, train_scores, test_scores = SVM_active_feature_selection_min_complexity(X_train, y_train, X_test, y_test, num_features=3, num_samples=num_samples, balance=balance)
    elif method == "min_cell":
        idx, num_samples_list, samples_global, train_errors, test_errors, train_scores, test_scores = SVM_active_feature_selection_min_cell(X_train, y_train, X_test, y_test, num_features=n, num_samples=num_samples)

    end = time.time()
    took = end - start

    # get names of selected indices
    selectedGenes = adata.var.index[idx].values
    selection = pd.DataFrame({"idx": idx, "selection": selectedGenes})

    return selection, took
