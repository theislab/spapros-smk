from selection_methods.gene_selection_shared import normalise
import numpy as np
import pandas as pd
import scanpy as sc
import sklearn as sk
import torch
from persist import PERSIST, ExpressionDataset, HurdleLoss


def preprocess_adata_persist(adata, ct_key="celltype", size_factors="size_factors"):
    """Preprocess adata for method PERSIST

    Arguments
    ----------
    adata: AnnData
        adata object with raw counts in adata.X and normalisation size factors in adata.obs, reduced to highly variable genes
    size_factors: str
        Size factors key
    """

    normalise(adata, key=size_factors)
    sc.pp.log1p(adata)
    # convert cell types to number ids
    adata.obs[ct_key] = pd.Categorical(adata.obs[ct_key]).codes 
    # save binarized data in a separate layer
    adata.layers['bin'] = (adata.X>0).astype(np.float32)
    
    return adata



def select_genes_persist(n, adata, ct_key="celltype", classification=False, train_size=0.8, max_nepochs=250, max_trials=100):
    """
    
    Parameter
    ---------
    classification: bool
        Wether to optimize for cell type classification. (if False we optimize for reconstruction, i.e. unsupervised case)
    train_size: float
        Ratio for training data size for train test split.
    max_nepochs: int
        Max number of epochs for elimination and selection.
    max_trials: int
        Number of trials for the elimination step before it throws an Error.
    
    
    """

    torch.manual_seed(777)
    np.random.seed(777)
    
    train_ind, val_ind = sk.model_selection.train_test_split(np.arange(adata.shape[0]), train_size=train_size, random_state=0)

    print("### PERSIST train test split:")
    print(f'\t {adata.shape[0]} total samples')
    print(f'\t {np.size(train_ind)} in training set')
    print(f'\t {np.size(val_ind)} in validation set')
    
    # These are views, so they do not take up memory
    adata_train = adata[train_ind,:]
    adata_val = adata[val_ind,:]
    
    # Initialize the dataset for PERSIST
    # Note: Here, data_train.layers['bin'] is a sparse array
    # data_train.layers['bin'].A converts it to a dense array
    if classification:
        train_dataset = ExpressionDataset(adata_train.layers['bin'].A, adata_train.obs[ct_key])
        val_dataset = ExpressionDataset(adata_val.layers['bin'].A, adata_val.obs[ct_key])
    else:
        train_dataset = ExpressionDataset(adata_train.layers['bin'].A, adata_train.X.A)
        val_dataset = ExpressionDataset(adata_val.layers['bin'].A, adata_val.X.A)
    
    
    # Use GPU device if available -- we highly recommend using a GPU!
    device = torch.device(torch.cuda.current_device() if torch.cuda.is_available() else 'cpu')
    print("### Device used (GPU highly recommended): ", device)
        
    # Set up the PERSIST selector
    selector = PERSIST(train_dataset,
                       val_dataset,
                       loss_fn=torch.nn.CrossEntropyLoss() if classification else HurdleLoss(),
                       device=device)


    # Ignore FutureWarning within selection block
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
    
        # Coarse removal of genes
        print('Starting initial elimination...')
        candidates, model = selector.eliminate(target=500, max_nepochs=max_nepochs, max_trials=max_trials)
        print('Completed initial elimination.')
        
        print('Selecting specific number of genes...')
        inds, model = selector.select(num_genes=n, max_nepochs=max_nepochs)
        persist_results = inds
        print('Done')
    
    genes = adata.var.iloc[inds].index.tolist()

    return genes
    