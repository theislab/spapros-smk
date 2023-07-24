import pandas as pd
import scanpy as sc
import time
import subprocess
import os

from selection_methods.gene_selection_shared import *


def preprocess_adata_selfe(adata, size_factors="size_factors", subsample=10000):
    """Preprocess adata for method selfE
    
    I actually just copied the preprocessing from scmer, we only included subsampling on top. Since SelfE always takes
    very long we use the same subsample size as for scmer.

    """

    normalise(adata, key=size_factors)
    sc.pp.log1p(adata)
    
    if subsample and (subsample < adata.n_obs):
        obs = np.random.choice(adata.n_obs, subsample, replace=False)
        adata = adata[obs]

    return adata


def select_genes_selfe(n, adata, tmp_dir, r_exe): #conda_env, env_dir, , size_factors="size_factors", subsample=10000
    #adata = preprocess_adata_selfe(adata, size_factors=size_factors, subsample=subsample)
    #adata = preprocess_adata(adata, size_factors=size_factors, subsample=10000)

    #print("Preprocess adata", flush=True)
    #preprocessed_adata_path = os.path.join(tmp_dir, "preprocessed_adata.h5ad")
    #print("Write adata to file", flush=True)
    #adata.write_h5ad(preprocessed_adata_path)
    #del adata  # need to delete the anndata because otherwise the h5ad file is locked and can't be read by the R script

    file_path = os.path.realpath(__file__)
    dir_path = os.path.dirname(file_path)
    r_script = dir_path + "/gene_selection_selfE.R"

    print("Start selection", flush=True)
    start = time.time()
    command = r_exe + " " + r_script + " " + \
              " ".join([str(n),  # 1
                        adata, #preprocessed_adata_path,  # 2
                        os.path.join(tmp_dir, "selfe_tmp.tsv"),  # 3
                        dir_path, # 4
                        #str(os.path.join(env_dir, conda_env))  # 4
                        ])
    print(command, flush=True)
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    #process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, bufsize=1)
    #for line in iter(process.stdout.readline, b''):
    #    print(line, flush=True)
    #process.stdout.close()    

    process.wait()
    end = time.time()
    print("Selection finished.", flush=True)
    took = end - start
    adata = sc.read(adata)
    idx = list(pd.read_csv(os.path.join(tmp_dir, "selfe_tmp.tsv"), sep="\t")['x'])
    selectedGenes = adata.var.index[idx]
    selection = pd.DataFrame({"idx": idx, "selection": selectedGenes})

    return selection, took
