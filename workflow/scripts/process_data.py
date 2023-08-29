import argparse
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse, csr_matrix
from util import normalise, get_highly_variable, subset_to_n_celltypes, sample_n_cells_per_ct, bootstrap_sample

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
    bootstrap_seed = int(params["bootstrap_seed"]) if params["bootstrap_seed"] != "None" else None
    
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
    
    # Subset to n cell types
    if n_cts is not None:
        adata = subset_to_n_celltypes(adata, n_cts, ct_key=ct_key, exact=False)
        
    # Subset each cell type to cells_per_ct
    if cells_per_ct is not None:
        adata = sample_n_cells_per_ct(adata, cells_per_ct, ct_key=ct_key, seed=cells_per_ct_seed)
        
    # Bootstrap sample of dataset, keep cell type proportions
    if bootstrap_seed is not None:
        adata = bootstrap_sample(adata, obs_key=ct_key, noise_level=1e-3, seed=bootstrap_seed, obs_names_unique=True)
    
    
    # Save data
    if not issparse(adata.X):
        adata.X = csr_matrix(adata.X)
    print(adata)
    adata.write(args.output)
    
                



if __name__ == '__main__':
    
    main()