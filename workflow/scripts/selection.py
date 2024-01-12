import os
import time
from pathlib import Path
import sys
import argparse
import pandas as pd
import scanpy as sc
from util import preprocess_adata, get_selection_df



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
    parser.add_argument('-s', '--save_specific_output', help='Save method specific output', required=False, type=str, default="False")
    
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
    save_specific_output = args.save_specific_output in ["True", True]
    
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
    selection, computation_time = run_selection(
        args.method, adata, n, ct_key, gene_key, proc, kwargs, 
        selection_csv, save_specific_output=save_specific_output
    )
    
    # Save selection
    df_selection = pd.DataFrame(index=adata.var.index, data={selection_id:False})
    df_selection[selection_id].iloc[selection.idx] = True
    df_info = pd.DataFrame(
        index=[selection_id],
        data={"time_seconds": [computation_time]}
    )

    df_selection.to_csv(selection_csv)
    df_info.to_csv(info_csv)
    
    
def run_selection(method, adata, n, ct_key, gene_key, proc, kwargs, selection_csv, save_specific_output=False):
    """
    """
    
    if gene_key is not None:
        adata = adata[:,adata.var[gene_key]]
        
    tmp_dir = Path(selection_csv.parent.parent , "tmp", selection_csv.stem)
    specific_dir = Path(selection_csv.parent.parent , "method_specific", selection_csv.stem)
    
    # SPAPROS
    if method == "spapros":
        import spapros as sp
        if proc:
            adata = preprocess_adata(adata)
        kwargs["n_min_markers"] = 0
        kwargs["genes_key"] = gene_key
        if save_specific_output:
            kwargs["save_dir"] = specific_dir
            
        selector = sp.se.ProbesetSelector(adata,ct_key,n=n,**kwargs,verbosity=0,n_jobs=-1)
        start = time.time()
        selector.select_probeset()
        computation_time = time.time() - start
        
        selection = get_selection_df(adata, selector.probeset.loc[selector.probeset["selection"]].index.tolist())
        
    # SPAPROS-CTO
    if method == "spaproscto":
        import spapros as sp
        if proc:
            adata = preprocess_adata(adata)
        kwargs["n_min_markers"] = 0
        kwargs["genes_key"] = gene_key
        if save_specific_output:
            kwargs["save_dir"] = specific_dir
            
        selector = sp.se.ProbesetSelector(adata,ct_key,n=n,n_pca_genes=0,**kwargs,verbosity=0,n_jobs=-1)
        start = time.time()
        selector.select_probeset()
        genes1 = selector.probeset.loc[selector.probeset["selection"]].index.tolist()
        selection1 = get_selection_df(adata, genes1)        
        
        # Run a 2nd selection if not enough genes were selected in the first round
        if len(genes1) < n:
            if save_specific_output:
                kwargs["save_dir"] = Path(specific_dir,"_run2")
            selector = sp.se.ProbesetSelector(
                adata[:,~adata.var_names.isin(genes1)],ct_key,n=n-len(genes1),n_pca_genes=0,**kwargs,verbosity=0,
                n_jobs=-1
            )
            selector.select_probeset()
            selection = get_selection_df(adata, genes1 + selector.probeset.loc[selector.probeset["selection"]].index.tolist())
        else:
            selection = selection1
            
        computation_time = time.time() - start
    
    # PCA
    elif method == "pca":
        import spapros as sp
        if proc:
            adata = preprocess_adata(adata)
            
        start = time.time()
        selection = sp.se.select_pca_genes(adata, n, inplace=False, verbosity=0, **kwargs)
        computation_time = time.time() - start
        
        selection = get_selection_df(adata, selection.loc[selection["selection"]].index.tolist())
    
    # DE
    elif method == "DE":
        import spapros as sp
        if proc:
            adata = preprocess_adata(adata)
            
        start = time.time()
        selection = sp.se.select_DE_genes(adata, n, obs_key=ct_key, inplace=False, verbosity=0, **kwargs)
        computation_time = time.time() - start
        
        selection = get_selection_df(adata, selection.loc[selection["selection"]].index.tolist())
    
    # SCGENEFIT
    elif method == "scgenefit":
        from selection_methods.gene_selection_scgenefit import select_genes_scgenefit, preprocess_adata_scgenefit
        if proc:
            adata = preprocess_adata_scgenefit(adata) #, **kwargs)
        kwargs["hierarchical"] = False
        kwargs["epsilon"] = 10
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
    elif method == "scpnmf":
        from selection_methods.gene_selection_scpnmf import select_genes_scpnmf, preprocess_adata_scpnmf
        if proc:
            adata = preprocess_adata(adata)
        if proc:
            adata = preprocess_adata_scpnmf(adata, subsample=10000)
        tmp_dir.mkdir(parents=True, exist_ok=True)
        #tmp_dir = os.path.join(out_dir, "tmp")
        #adata = os.path.join(general_params["data_path"],general_params["dataset"])  # need to delete the anndata because otherwise the h5ad file is locked and can't be read by the R script
        #if not os.path.exists(tmp_dir):
        #    os.umask(0)
        #    os.makedirs(tmp_dir, 0o777)
        #    os.chmod(tmp_dir, 0o777)
        adata_path = os.path.join(tmp_dir, "adata_tmp.h5ad")
        adata.write(adata_path)
        r_exe = "Rscript" #"/usr/bin/Rscript"
        selection, computation_time = select_genes_scpnmf(
            n, adata_path, r_exe, output_path=os.path.join(tmp_dir, "scpnmf_tmp.tsv"), **kwargs
        )

    # SCMER
    elif method == "scmer":
        from selection_methods.gene_selection_scmer import select_genes_scmer, preprocess_adata_scmer
        if proc:
            adata = preprocess_adata_scmer(adata, n_pcs=30, subsample=10000)
        kwargs["n_threads"] = -1
        selection_scmer, computation_time = select_genes_scmer(n, adata, **kwargs)
        
        if len(selection_scmer) == n:
            selection = selection_scmer
        else:
            print("SCMER selected ", len(selection_scmer), "genes, but ", n, "were requested. Adding/removing based on pca scores.")
            from util import select_pca_genes
            pca_scores = select_pca_genes(adata, 0, inplace=False)
            
            pca_scores['idx'] = range(adata.shape[1])
            pca_scores["scmer_selection"] = False
            pca_scores.loc[selection_scmer["selection"].values, "selection"] = True
            pca_scores.loc[selection_scmer["selection"].values, "scmer_selection"] = True
            pca_scores = pca_scores.sort_values(["selection","selection_score"],ascending=[False,False])
            # set first n values of "selection" to True and rest to False
            pca_scores["selection"] = False
            pca_scores.loc[pca_scores.index[:n], "selection"] = True
            
            if save_specific_output:
                Path(specific_dir).mkdir(parents=True, exist_ok=True)
                pca_scores.to_csv(Path(specific_dir) / "selection_adjusted_to_n_genes.csv")
            
            selection = get_selection_df(adata, pca_scores.loc[pca_scores["selection"]].index.tolist())

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
    
    # selfE
    elif method == "selfe":
        from selection_methods.gene_selection_selfe import select_genes_selfe, preprocess_adata_selfe
        if proc:
            adata = preprocess_adata_selfe(adata, subsample=10000)
        tmp_dir.mkdir(parents=True, exist_ok=True)
        #conda_env = config["venv"]
        #if not os.path.exists(tmp_dir):
        #    os.umask(0)
        #    os.makedirs(tmp_dir, 0o777)
        #    os.chmod(tmp_dir, 0o777)
        adata_path = os.path.join(tmp_dir, "adata_tmp.h5ad")
        adata.write(adata_path)
        r_exe = "Rscript" #"/usr/bin/Rscript"        
        selection, computation_time = select_genes_selfe(n, adata_path, tmp_dir, r_exe)

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

