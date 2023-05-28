import os
from pathlib import Path
import sys
import argparse
import pandas as pd
import scanpy as sc
from util import preprocess_adata



# Get arguments
def get_args():
    """Get input arguments
    """
    parser = argparse.ArgumentParser(description="Run selection")
    
    parser.add_argument('-d', '--data', help='Input h5ad file name', required=True, type=str)
    parser.add_argument('-o', '--output', help='Output directory', required=True, type=str)
    parser.add_argument('-m', '--method', help='Selection method', required=True, type=str)
    parser.add_argument('-id', '--id', help='Selection id', required=True, type=int)
    parser.add_argument('-f', '--id_file', help='Csv file with ids and parameter sets', required=True, type=str)
    
    return parser.parse_args()
    
    
def main():
    """
    """
    
    # Get input arguments
    args = get_args()
    
    # Define selection id and output files
    selection_id = args.method + "_" + str(args.id)
    selection_csv = Path(args.output)
    info_csv = Path(args.output.replace(".csv","_info.csv"))
    
    # Get parameters
    params = pd.read_csv(args.id_file, index_col=0)
    for col in params.columns:
        params.loc[params[col].isnull(), col] = "None"
    params = params.loc[args.id].to_dict()
    
    n = int(params["n"]) if params["n"] != "None" else None
    ct_key = str(params["ct_key"]) if params["ct_key"] != "None" else None
    gene_key = str(params["gene_key"]) if params["gene_key"] != "None" else None
    proc = True if (params["method_specific_processing"] in ["True", True]) else False
    
    kwargs = {} # we might want to extend the pipeline to method specific parameters
    
    # Load data
    adata = sc.read(args.data)
    if proc and (adata.uns["processing"] != "raw"):
        raise ValueError("Data needs to be raw if method specific processing is enabled, in the config set either `processing` to None or `method_specific_processing` to False")
    
    # Reduce to pre filtered genes
    if gene_key:
        adata = adata[:,adata.var[gene_key]]
    
    # Run selection
    selection, computation_time = run_selection(args.method, adata, n, ct_key, gene_key, proc, kwargs, selection_csv)
    
    # Save selection
    df_selection = pd.DataFrame(index=adata.var.index, data={selection_id:False})
    df_selection[selection_id].iloc[selection.idx] = True
    df_info = pd.DataFrame(
        index=[selection_id],
        data={"time_seconds": [computation_time]}
    )

    df_selection.to_csv(selection_csv)
    df_info.to_csv(info_csv)
    
    
def run_selection(method, adata, n, ct_key, gene_key, proc, kwargs, selection_csv):
    """
    """
    
    adata = adata[:,gene_key]
    
    # SPAPROS
    if method == "spapros":
        import spapros as sp
        if proc:
            adata = preprocess_adata(adata)
        kwargs["n_min_markers"] = 0
        kwargs["genes_key"] = gene_key
        selector = sp.se.ProbesetSelector(adata,ct_key,n=n,**kwargs,verbosity=0,n_jobs=-1)
        selector.select_probeset()
        selection = selector.probeset
    
    # PCA
    elif method == "pca":
        import spapros as sp
        if proc:
            adata = preprocess_adata(adata)        
        selection = sp.se.select_pca_genes(adata, n, inplace=False, verbosity=0, **kwargs)
    
    # DE
    elif method == "DE":
        import spapros as sp
        if proc:
            adata = preprocess_adata(adata)        
        selection = sp.se.select_DE_genes(adata, n, obs_key=ct_key, inplace=False, verbosity=0, **kwargs)
    
    # SCGENEFIT
    elif method == "scgenefit":
        from selection_methods.gene_selection_scgenefit import select_genes_scgenefit, preprocess_adata_scgenefit
        if proc:
            adata = preprocess_adata_scgenefit(adata) #, **kwargs)
        kwargs["hierarchical"] = True
        kwargs["epsilon"] = 0.5          
        selection, computation_time = select_genes_scgenefit(n,adata=adata,label=ct_key,**kwargs)

    # NSFOREST
    elif method == "nsforest":
        from selection_methods.gene_selection_nsforest import gene_selection_nsforest
        if proc:
            print("############## SPECIFIC PREPROCESSING ##############")            
            preprocess_adata(adata)
        # labl = kwargs["label"]
        # n_clust = adata.obs[labl].cat.categories.shape[0]
        selection, computation_time = gene_selection_nsforest(n, adata=adata, label=ct_key, **kwargs)  # int(np.ceil(n / n_clust))

    ## SCPNMF
    #elif method == "scpnmf":
    #    from selection_methods.gene_selection_scpnmf import select_genes_scpnmf
    #    tmp_dir = os.path.join(out_dir, "tmp")
    #    conda_env = config["venv"]
    #    adata = os.path.join(general_params["data_path"],general_params["dataset"])  # need to delete the anndata because otherwise the h5ad file is locked and can't be read by the R script
    #    if not os.path.exists(tmp_dir):
    #        os.umask(0)
    #        os.makedirs(tmp_dir, 0o777)
    #        os.chmod(tmp_dir, 0o777)
    #    # r_exe = "/usr/bin/Rscript"
    #    selection, computation_time = select_genes_scpnmf(n, adata, output_path=os.path.join(tmp_dir, "scpnmf_tmp.tsv"), **kwargs, conda_env=conda_env)
    #    adata = sc.read(adata)

    # SCMER
    elif method == "scmer":
        from selection_methods.gene_selection_scmer import select_genes_scmer, preprocess_adata_scmer
        if proc:
            adata = preprocess_adata_scmer(adata, n_pcs=30, subsample=10000)
        kwargs["n_threads"] = -1
        selection, computation_time = select_genes_scmer(n, adata, **kwargs)

    # SMASH
    elif method == "smash":
        from selection_methods.gene_selection_smash import select_genes_smash, preprocess_adata_smash
        if proc:
            adata = preprocess_adata_smash(adata.copy(), ct_key)
        kwargs["method"] = "XGB"
        selection, computation_time = select_genes_smash(n, adata, ct_key, **kwargs)
        #if 'method' in kwargs:
        #    name = name + "_" + kwargs["method"]
        #else:
        #    name = name + "_DNN"

    # ASFs
    elif method == "asfs":
        from selection_methods.gene_selection_asfs import select_genes_asfs
        if proc:
            adata = preprocess_adata(adata)
        kwargs["num_samples"] = 300
        selection, computation_time = select_genes_asfs(n, adata, ct_key, **kwargs)

    # geneBasis
    elif method == "genebasis":
        from selection_methods.gene_selection_genebasis import select_genes_genebasis
        if proc:
            adata = preprocess_adata(adata)
        tmp_dir = Path(selection_csv.parent.parent , "tmp", selection_csv.stem)
        tmp_dir.mkdir(parents=True, exist_ok=True)
        #tmp_dir = os.path.join(out_dir, "tmp")
        #conda_env = config["venv"]
        #if not os.path.exists(tmp_dir):
        #    os.umask(0)
        #    os.makedirs(tmp_dir, 0o777)
        #    os.chmod(tmp_dir, 0o777)
        selection, computation_time = select_genes_genebasis(n=n,
                                                             adata=adata,
                                                             tmp_dir=tmp_dir,
                                                             #conda_env=conda_env,
                                                             **kwargs)
    
    ## selfE
    #elif method == "selfe":
    #    from selection_methods.gene_selection_selfe import select_genes_selfe
    #    tmp_dir = os.path.join(out_dir, "tmp")
    #    conda_env = config["venv"]
    #    if not os.path.exists(tmp_dir):
    #        os.umask(0)
    #        os.makedirs(tmp_dir, 0o777)
    #        os.chmod(tmp_dir, 0o777)
    #    selection, computation_time = select_genes_selfe(n=n,
    #                                                     adata=adata,
    #                                                     tmp_dir=tmp_dir,
    #                                                     conda_env=conda_env,
    #                                                     **kwargs)

    # COSG
    elif method == "cosg":
        from selection_methods.gene_selection_cosg import select_genes_cosg
        if proc:
            adata = preprocess_adata(adata)
        kwargs["groupby"] = ct_key
        selection, computation_time = select_genes_cosg(n, adata, **kwargs)

    # Triku
    elif method == "triku":
        from selection_methods.gene_selection_triku import select_genes_triku, preprocess_adata_triku
        if proc:
            adata = preprocess_adata_triku(adata)
        selection, computation_time = select_genes_triku(n, adata, **kwargs)
    
    
    return selection, computation_time
    
    


if __name__ == "__main__":

    main()

