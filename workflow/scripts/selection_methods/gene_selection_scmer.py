from selection_methods.gene_selection_shared import *
import scmer
import scanpy as sc
import numpy as np
import pandas as pd
import time



def preprocess_adata_scmer(adata, size_factors="size_factors", n_pcs=100, subsample=10000):
    """Preprocess adata for method SCMER

    Arguments
    ----------
    adata: AnnData
        adata object with raw counts in adata.X and normalisation size factors in adata.obs, reduced to highly variable genes
    label: str
        Celltype clustering key
    n_pcs: Number of principal components to compute.

    the recommended preprocessing workflow for scmer is:

    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.normalize_total(adata, target_sum=1e4)
    # sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    sc.tl.pca(adata, svd_solver='arpack', n_comps=100)
    sc.pp.neighbors(adata, n_pcs=100)

    But for camparability, we use size factors like with the other methods

    """


    normalise(adata, key=size_factors)
    sc.pp.log1p(adata)
    
    np.random.seed(0)
    
    if subsample and (subsample < adata.n_obs):
        obs = np.random.choice(adata.n_obs, subsample, replace=False)
        adata = adata[obs]

    return adata


def select_genes_scmer(n, adata, n_pcs=100, size_factors="size_factors", w='ones', lasso=0.0001, perplexity=30.0, use_beta_in_Q=True,
                       max_outer_iter=5, max_inner_iter=20, owlqn_history_size=32, eps=1e-12, verbosity=2,
                       torch_precision=32, torch_cdist_compute_mode='use_mm_for_euclid_dist', t_distr=True, n_threads=1,
                       use_gpu=False, pca_seed=0, ridge=0.0, _keep_fitting_info=False, X_teacher=None, batches=None,
                       P=None, beta=None, must_keep=None, min_lasso=1e-8, max_lasso=0.01, tolerance=0,
                       smallest_log10_fold_change=0.1, max_iter=100, n_threads_fit=6, return_P_beta=False, subsample=10000):
    """Select genes via SCMER

    Arguments
    -----------
    n: int
        Number of genes to select
    adata: AnnData
        adata object with normalized and logarithmized counts in adata.X
    w: Union[float, str, list, numpy.ndarray] (default: 'ones')
        initial value of w, weight of each marker. Acceptable values are ‘ones’ (all 1), ‘uniform’ (random [0, 1] values), float numbers (all set to that number), or a list or numpy array with specific numbers.
    lasso: float (default: 0.0001)
        lasso strength (i.e., strength of L1 regularization in elastic net)
    n_pcs: int (default: None)
        Number of PCs used to generate P matrix. Skip PCA if set to None.
    size_factors: (default: "size_factors")
        size factors for normalizing
    perplexity: float (default: 30.0)
        perplexity of t-SNE modeling
    use_beta_in_Q: bool (default: True)
        whether to use the cell specific sigma^2 calculated from P in Q. (1 / beta)
    max_outer_iter: int (default: 5)
        number of iterations of OWL-QN
    max_inner_iter: int (default: 50)
        number of iterations inside OWL-QN
    owlqn_history_size: int (default: 100)
        history size for OWL-QN.
    eps: float (default: 1e-12)
        epsilon for considering a value to be 0.
    verbosity: int (default: 2)
        verbosity level (0 ~ 2).
    torch_precision: Union[int, str, torch.dtype] (default: 32)
        The dtype used inside torch model. By default, tf.float32 (a.k.a. tf.float) is used. However, if precision become an issue, tf.float64 may be worth trying. You can input 32, “32”, 64, or “64”.
    torch_cdist_compute_mode: str (default: 'use_mm_for_euclid_dist')
        cdist_compute_mode: compute mode for torch.cdist. By default, “use_mm_for_euclid_dist” to (daramatically) improve performance. However, if numerical stability became an issue, “donot_use_mm_for_euclid_dist” may be used instead. This option does not affect distances computed outside of pytorch, e.g., matrix P. Only matrix Q is affect.
    t_distr: bool (default: True)
        By default, use t-distribution (1. / (1. + pdist2)) for Q. Use Normal distribution instead (exp(-pdist2)) if set to False. The latter one is not stable.
    n_threads: int (default: 1)
        number of threads (currently only for calculating P and beta)
    use_gpu: bool (defualt: False)
        whether to use GPU to train the model.
    pca_seed: int (default: 0)
        random seed used by PCA (if applicable)
    ridge: float (default: 0.0)
        ridge strength (i.e., strength of L2 regularization in elastic net)
    _keep_fitting_info: bool (default: False)
        if True, write similarity matrix P to self.P and PyTorch model to self.model
    X_teacher: matrix (default: None)
        get target similarities from this dataset
    batches: list (default: None)
        (optional) batch labels
    P: matrix (default: None)
        The P matrix, if calculated in advance
    beta: (default: None)
        The beta associated with P, if calculated in advance
    must_keep: vector (default: None)
        A bool vector indicating if a feature must be kept. Those features will have a fixed weight 1.
    min_lasso: float (default: 1e-8)
        nowhere mentioned
    max_lasso: float (default: 0.01)
        nowhere mentioned
    tolerance: int (default: 0)
        nowhere mentioned
    smallest_log10_fold_change: float (default: 0.1)
        nowhere mentioned
    max_iter: int (default: 100)
        nowhere mentioned
    n_threads_fit: int
        number of threads (currently only for calculating P and beta) during fitting
    return_P_beta: bool (default: False)
        controls what to return
        :param size_factors:
    subsample: int (default: 10000)
        Take random subsample of dataset. Scmer computation time grows quadratic therefore the authors recommend to 
        subsample the dataset to 5000 - 10000 cells.
    """

    #adata = preprocess_adata_scmer(adata, size_factors, subsample=subsample)

    # np.ndarray needed
    adata_cp = adata.copy()
    if not type(adata.X) == np.ndarray:
        adata_cp.X = adata.copy().X.toarray()

    start = time.time()

    # specify model
    model = scmer.UmapL1(w=w, lasso=lasso, n_pcs=n_pcs, perplexity=perplexity, use_beta_in_Q=use_beta_in_Q,
                         max_outer_iter=max_outer_iter, max_inner_iter=max_inner_iter,
                         owlqn_history_size=owlqn_history_size, eps=eps, verbosity=verbosity,
                         torch_precision=torch_precision, torch_cdist_compute_mode=torch_cdist_compute_mode,
                         t_distr=t_distr, n_threads=n_threads, use_gpu=use_gpu, pca_seed=pca_seed,
                         ridge=ridge)  # , _keep_fitting_info=_keep_fitting_info)

    # like model.fit() but automatically find proper lasso strength that returns the preferred number of markers

    model = model.tune(n, X=adata_cp.X, X_teacher=X_teacher, batches=batches, P=P, beta=beta, must_keep=must_keep,
                       perplexity=perplexity, n_pcs=n_pcs, w=w, min_lasso=min_lasso, max_lasso=max_lasso, tolerance=tolerance,
                       smallest_log10_fold_change=smallest_log10_fold_change, max_iter=max_iter, return_P_beta=return_P_beta,
                       n_threads=n_threads_fit)

    end = time.time()
    took = end - start

    selectedGenes = adata_cp.var_names[model.w > 0]

    # get index in adata for selected genes
    adata.var['idx'] = range(adata.shape[1])
    idx = [int(adata.var['idx'][adata.var.index == x].values[0]) for x in selectedGenes]
    selection = pd.DataFrame({"idx": idx, "selection": selectedGenes})

    return selection, took

