from selection_methods.NS_Forest_3 import *
from selection_methods.gene_selection_shared import *
import pandas as pd
import time



def gene_selection_nsforest(n, adata, label, rfTrees=1000, threads=-1, Median_Expression_Level=0, Genes_to_testing=3, betaValue=0.5):
    
    #adata = preprocess_adata(adata)
    
    # nsforest doesn't work with special characters, use simple index for genes instead
    initial_var_names = adata.var_names.tolist().copy()
    int_to_str = lambda x: ''.join([chr(int(i)+65) for i in str(x)])
    dummy_index = [int_to_str(i) for i in range(adata.n_vars)]
    adata.var_names = dummy_index
    mapping = {idx:gene for idx, gene in zip(dummy_index,initial_var_names)}

    start = time.time()
    binary_table = coreAnalysis(n, adata, clusterLabelcolumnHeader=label, rfTrees=rfTrees, threads=threads,
                                Median_Expression_Level=Median_Expression_Level, Genes_to_testing=Genes_to_testing, betaValue=betaValue)
    # binary_table = NS_Forest(adata, clusterLabelcolumnHeader=label, rfTrees=rfTrees, Median_Expression_Level=Median_Expression_Level, Genes_to_testing=n, betaValue=betaValue)

    end = time.time()
    took = end - start
    
    # Map simple indices back to gene names
    adata.var_names = initial_var_names
    binary_table = binary_table.rename(index=mapping)

    # b = binary_table[::-n_per_clust]
    # n_clust = adata.obs[label].cat.categories.shape[0]
    # diff = (n_clust * n_per_clust) - n
    # to_remove = b[:diff].index.to_list()
    # binary_table = binary_table.drop(to_remove, axis=0)

    # selectedGenes = binary_table.index.tolist()

    selectedGenes = []
    clusters = binary_table.clusterName.unique()
    remaining_clusters = binary_table.clusterName.unique()
    while len(selectedGenes) < n and len(remaining_clusters) > 0:
        for cluster in clusters:
            genes_of_cluster = binary_table[binary_table.clusterName == cluster].sort_values(by="0_informationGain", ascending=False)
            if len(genes_of_cluster) == 0:
                continue
            best_gene_of_cluster = genes_of_cluster.iloc[0, :].name
            if best_gene_of_cluster not in selectedGenes and len(selectedGenes) < n:
                selectedGenes.append(best_gene_of_cluster)
            binary_table = binary_table.drop(best_gene_of_cluster, axis=0)
            remaining_clusters = binary_table.clusterName.unique()

    adata.var['idx'] = range(adata.shape[1])
    #idx = [int(adata.var['idx'][adata.var.gene_ids == x].values[0]) for x in selectedGenes]
    idx = [int(adata.var['idx'][adata.var_names == x].values[0]) for x in selectedGenes]    
    names = adata.var.index.values[idx]
    selection = pd.DataFrame({"idx": idx, "selection": names})

    # selected = pd.Series([x in selectedGenes for x in adata.var.gene_ids])
    # selected.index = adata.var.index
    # selection = pd.DataFrame(selected, columns=["selection"])

    return selection, took
