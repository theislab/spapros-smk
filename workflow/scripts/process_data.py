import argparse
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse, csr_matrix, sparsefuncs

def get_args():
    """Get input arguments"""
    parser = argparse.ArgumentParser(description='Process dataset')
    
    parser.add_argument('-d', '--data', help='Input h5ad file name', required=True, type=str)
    parser.add_argument('-o', '--output', help='Output h5ad file name', required=True, type=str)
    parser.add_argument('-id', '--id', help='Dataset processing id', required=True, type=int)
    parser.add_argument('-f', '--id_file', help='Csv file with ids and parameter sets', required=True, type=str)
    
    return parser.parse_args()

def main():
    """Run data processing"""
    args = get_args()
    print(args.data)
    print(args.output)
    print(args.id)
    print(args.id_file)
    
    # Get parameters
    params = pd.read_csv(args.id_file, index_col=0)
    for col in params.columns:
        params.loc[params[col].isnull(), col] = "None"
    params = params.loc[args.id].to_dict()
    
    processing = str(params["processing"]) if params["processing"] != "None" else None
    ct_key = params["ct_key"] if params["ct_key"] != "None" else None
    n_cts = int(params["n_cts"]) if params["n_cts"] != "None" else None
    cells_per_ct_seed = int(params["cells_per_ct_seed"]) if params["cells_per_ct_seed"] != "None" else None
    cells_per_ct = int(params["cells_per_ct"]) if params["cells_per_ct"] != "None" else None
    
    # Load data
    adata = sc.read(args.data)
    print(adata)
    
    # Process data
    # Save size_factors based on normalize total (needed in case method specific normalization is done in the selection)
    if "size_factors" not in adata.obs.columns:
        a = adata.copy()
        tmp = sc.pp.normalize_total(a,target_sum=1e6,inplace=False)
        adata.obs["size_factors"] = tmp['norm_factor']/10**6
        
    if "highly_variable" not in adata.var.columns:
        get_highly_variable(adata,N=8000,key_added='highly_variable')
    
    if processing == "lognorm":
        normalise(adata)
        sc.pp.log1p(adata)
        adata.uns["processing"] = "lognorm"
    else:
        adata.uns["processing"] = "raw"
    
    # Save data
    if not issparse(adata.X):
        adata.X = csr_matrix(adata.X)
    print(adata)
    adata.write(args.output)
    


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

def get_highly_variable(adata,N=8000,key_added='highly_variable'):
    """Compute N highly variable genes of given dataset for given normalisation
    
    A copy of the raw data is normalised and log1p transformed to extract N highly variable genes.
    """
    a = adata.copy()
    a.X /= a.obs['size_factors'].values[:,None]
    sc.pp.log1p(a)
    sc.pp.highly_variable_genes(a, n_top_genes=N, flavor='cell_ranger', inplace=True)
    adata.var[key_added] = a.var['highly_variable']


if __name__ == '__main__':
    
    main()