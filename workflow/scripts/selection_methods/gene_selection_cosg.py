from selection_methods.gene_selection_shared import *
import time
import cosg
import pandas as pd


def select_genes_cosg(n, adata, size_factors="size_factors", **kwargs):
    """
    Parameters
    ----------
    n: int
        Number of genes to select
    adata: AnnData
        adata object with normalized and logarithmized counts in adata.X
    size_factors: String (default: "size_factors")
        size_factors for normalization
    kwargs:
        groupby: String (default: 'CellTypes')
            The key of the cell groups in .obs, the default value is set to 'CellTypes'.
        groups: Union[typing_extensions.Literal['all'], Iterable[str]] (defaut: 'all')
            Subset of cell groups, e.g. [`'g1'`, `'g2'`, `'g3'`], to which comparison shall be restricted. The default value is 'all', and all groups will be compared.
        mu: int (default: 1)
            The penalty restricting marker genes expressing in non-target cell groups. Larger value represents more strict restrictions. mu should be >= 0, and by default, mu = 1.
        remove_lowly_expressed: Boolean (default: False)
            If True, genes that express a percentage of target cells smaller than a specific value (`expressed_pct`) are not considered as marker genes for the target cells. The default value is False.
        expressed_pct: Union[float, NoneType] (default: 0.1)
            When `remove_lowly_expressed` is set to True, genes that express a percentage of target cells smaller than a specific value (`expressed_pct`) are not considered as marker genes for the target cells. The default value for `expressed_pct`
            is 0.1 (10%).
        key_added: Union[str, NoneType] (default: None)
            The key in `adata.uns` information is saved to.
        use_raw: Boolean (default: True)
            Use `raw` attribute of `adata` if present.
        layer: Union[str, NoneType] (default: None)
            Key from `adata.layers` whose value will be used to perform tests on.
        reference: String (default: 'rest')
            If `'rest'`, compare each group to the union of the rest of the group.
            If a group identifier, compare with respect to this group
        copy: Boolean (default: False)
    """

    #adata = preprocess_adata(adata, size_factors=size_factors)

    start = time.time()

    new_adata = cosg.cosg(adata, n_genes_user=n, copy=True, **kwargs)
    selectedGenes = list(set(sum([list(x) for x in new_adata.uns['rank_genes_groups']['names']], [])))[:n]

    end = time.time()
    took = end - start

    # get names of selected indices
    adata.var['idx'] = range(adata.shape[1])
    idx = [int(adata.var['idx'][adata.var.index == x]) for x in selectedGenes]
    selection = pd.DataFrame({"idx": idx, "selection": selectedGenes})

    return selection, took
