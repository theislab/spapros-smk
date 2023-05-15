import argparse
import pandas as pd
import scanpy as sc

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
    
    ct_key = params["ct_key"] if params["ct_key"] != "None" else None
    n_cts = int(params["n_cts"]) if params["n_cts"] != "None" else None
    cells_per_ct_n_seeds = int(params["cells_per_ct_n_seeds"]) if params["cells_per_ct_n_seeds"] != "None" else None
    cells_per_ct = int(params["cells_per_ct"]) if params["cells_per_ct"] != "None" else None
    
    # Load data
    adata = sc.read(args.data)
    print(adata)
    # Process data
    
    # Save data
    print(adata)
    adata.write(args.output)
    



if __name__ == '__main__':
    
    main()