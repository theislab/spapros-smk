import smashpy
sm = smashpy.smashpy()
from scipy.sparse import issparse
from selection_methods.gene_selection_shared import *
import scanpy as sc
import numpy as np
import pandas as pd
import time


def preprocess_adata_smash(adata, label):
    """Preprocess adata for method smash

    Arguments
    ---------
    adata: AnnData
        adata object with raw counts in adata.X, reduced to highly variable genes
    label: str
        Celltype clustering key
    """

    if issparse(adata.X):

        # norm, log, scale:
        adata.layers["counts"] = adata.X.copy().toarray()

        # normalise and save the data into the layer
        sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)  # deprecated, use normilize_total
        adata.layers["norm"] = adata.X.copy().toarray()

        # logarithmise and save the data into the layer
        sc.pp.log1p(adata)
        adata.layers["log"] = adata.X.copy().toarray()

        # save in adata.raw.X normilise and logarithm data
        adata.raw = adata.copy()

        # scale adata, attention this makes adata.X dense !
        sc.pp.scale(adata, max_value=10)
        adata.layers["scale"] = adata.X.copy()

    else:

        # norm, log, scale:
        adata.layers["counts"] = adata.X.copy()

        # normalise and save the data into the layer
        sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)  # deprecated, use normilize_total
        adata.layers["norm"] = adata.X.copy()

        # logarithmise and save the data into the layer
        sc.pp.log1p(adata)
        adata.layers["log"] = adata.X.copy()

        # save in adata.raw.X normalize and logarithmize data
        adata.raw = adata.copy()
        sc.pp.scale(adata, max_value=10)  # attention this makes adata.X dense !
        adata.layers["scale"] = adata.X.copy()

    # removing general genes
    adata = sm.remove_general_genes(adata)

    # removing genes expressed in less than 30% within groups (cell types)
    adata = sm.remove_features_pct(adata, group_by=label, pct=0.3)

    # rmoving genes expressed in more than 50% in a given group where genes are expressed for more 75% within a given group
    adata = sm.remove_features_pct_2groups(adata, group_by=label, pct1=0.75, pct2=0.5)

    # inverse pca
    adata = sm.scale_filter_features(adata, n_components=None, filter_expression=True)

    return adata


def select_genes_smash(n, adata, label, method="DNN", model=None, test_size=0.2, balance=True, verbose=False,
                       save=False, pct=0.05):

    """Select genes via SMaSH

    Arguments
    ----------
    n: int
        Number of genes to select
    adata: AnnData
        adata object with normalized and logarithmized counts in adata.X
    label: String
        name of adata.obs where a list of labels can be found: eg. cell type, cell state or
        list with labels (N labels, one per point (gene names))
    method: String (default 'DNN')
        'DNN', 'RF', 'BRF' or
    model : Keras model (default: None)
    test_size : float (default: 0.2)
    balance : bool (default: True)
    verbose : bool (default: True)
    save : bool (default: True)
    pct: float (default: 0.05)
        which fraction of test and train sets to use, if None, use full 0.2 for test and 0.8 for training
    """

    #adata_cp = preprocess_adata_smash(adata.copy(), label)
    adata_cp = adata.copy()

    start = time.time()

    # feedforward deep feedforward neural network (DNN)
    if method == "DNN":
        # Applying DNN to adata.X
        sm.DNN(adata_cp, group_by=label, model=model, test_size=test_size, balance=balance, verbose=verbose, save=save)
        # shapley value ranking
        selectedGenes, selectedGenes_dict = sm.run_shap(adata_cp, group_by=label, model=None, verbose=False, pct=pct,
                                                        restrict_top=('global', n))

    else:
        # random forest (RF)
        if method == "RF":
            # Applying ensemble_learning to adata.X
            clf = sm.ensemble_learning(adata_cp, group_by=label, classifier="RandomForest", test_size=test_size,
                                       balance=balance, verbose=verbose, save=save)

        # balanced random forest (BRF)
        elif method == 'BRF':
            # Applying ensemble_learning to adata.X
            clf = sm.ensemble_learning(adata_cp, group_by=label, classifier="BalancedRandomForest", test_size=test_size,
                                       balance=balance, verbose=verbose, save=save)

        # XGBoost (XGB)
        elif method == 'XGB':
            # Applying ensemble_learning to adata.X
            clf = sm.ensemble_learning(adata_cp, group_by=label, classifier="XGBoost", test_size=test_size,
                                       balance=balance, verbose=verbose, save=save)

        # select top n genes according to gini importance
        selectedGenes, selectedGenes_dict = sm.gini_importance(adata_cp, clf, group_by=label, verbose=verbose,
                                                               restrict_top=('global', n))

    end = time.time()
    took = end - start

    # get index in adata for selected genes
    adata.var['idx'] = range(adata.shape[1])
    idx = [int(adata.var['idx'][adata.var.index == x]) for x in selectedGenes]
    selection = pd.DataFrame({"idx": idx, "selection": selectedGenes})

    return selection, took
