import numpy as np
import pandas as pd

from pyroe import load_fry
import scanpy as sc
import decoupler as dc
import scvelo as scv

import umap

import matplotlib.pyplot as plt
import sys


sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor="white")


def fill_unmapped(name_df):
    '''Some ENTREZ IDs are not mapped to gene symbols. Fill them in or else downstream stuff complains'''
    na_i = name_df["GENE_SYMBOL"].isna()
    filler = [f"UNDEF{x}" for x in range(na_i.sum())]

    name_df.loc[name_df["GENE_SYMBOL"].isna(), "GENE_SYMBOL"] = filler

    return name_df


def number_duplicates(arr):
    '''Number duplicate entries so things don't complain'''
    cat_arr = pd.Series([""]*len(arr))
    dup_i = arr.duplicated()

    cat_arr[dup_i] = ".2"
    return arr.str.cat(cat_arr)


if __name__ == "__main__":
    # load count matrix
    srr_id = str(sys.argv[1])
    adata = load_fry(f"/nobackup1/ericjlu/quantsv3/{srr_id}_quant_v3/af_quant/", output_format = "velocity")
    
    # Load gene name mapping
    name_df = pd.read_csv("./gene_name_mapping.csv")
    name_df = fill_unmapped(name_df)
    
    
    # Basic filtering
    initial_cts = adata.shape[0]
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    
    # Annotate the group of mitochondrial genes as 'mt'
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    # Filter cells following standard QC criteria.
    adata = adata[adata.obs.n_genes_by_counts < 2500, :]
    adata = adata[adata.obs.pct_counts_mt < 5, :]
    
    # Normalize the data
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    adata.layers['log_norm'] = adata.X.copy()
    
    print(f"Pass rate: {np.around(adata.shape[0]/initial_cts*100, 2)}%")
    
    
    # Finagle with the namings
    merged = adata.var.merge(name_df, how="inner", left_index=True, right_on="ENSG_ID")
    merged.index = merged["GENE_SYMBOL"]
    merged.index = number_duplicates(pd.Series(merged.index))
    merged.index.name = "gene_ids"
    merged = merged.drop(["ENSG_ID", "GENE_SYMBOL"], axis=1)
    
    adata.var = merged
    adata.var.index = number_duplicates(pd.Series(adata.var.index))
    
    sc.pp.highly_variable_genes(adata)
    sc.pp.neighbors(adata)
    scv.tl.umap(adata, random_state=1)
    umap_x = adata.obsm["X_umap"]
    
    # Full dynamics
    scv.tl.recover_dynamics(adata, n_jobs=1, var_names="highly_variable")
    # Expedited
    scv.tl.velocity(adata, n_jobs=1)
    
    
    # These genes come from the nature paper
    ct_markers = {
        "Mesenchymal cells": ["DCN", "COL11A1", "PDGFRA","PDGFRB"],
        "Epithelial cells": ["EPCAM", "KRT8", "KRT10", "KRT18", "KRT19"],
        "Smooth muscle cells": ["ACTA2"],
        "Erythrocytes": ["HBB", "GYPA"],
        "Endothelial cells": ["CLDN5", "PECAM1", "CD34", "ESAM"],
        "Mast cells": ["TPSB2"],
        "Myeloid cells": ["LYZ", "CD14","C1QA", "CLEC10A"],
        "T/NK cells": ["CD2", "CD3D", "CD3E", "CD3G", 'CD8A', "CCL5"],
        "B/Plasma cells": ["JCHAIN", "CD79A"]
    }
    
    cts_df = adata.to_df()
    mgs = np.concatenate(list(ct_markers.values()))
    avail_mgs = np.intersect1d(cts_df.columns, mgs)
    print("Markers filtered out:",  ', '.join(np.setdiff1d(mgs, avail_mgs)))
    
    # Record filtered-out genes
    mg_df = cts_df[avail_mgs]
    filtered_mgs = list(np.setdiff1d(mgs, avail_mgs))
    
    # Remove filtered-out genes from ct_markers
    for k, v in ct_markers.items():
        for i in range(len(filtered_mgs)):
            if filtered_mgs[i] in ct_markers[k]:
                ct_markers[k].remove(filtered_mgs[i])
    
    # Louvain clustering
    sc.tl.louvain(adata)
    
    # Ranking and testing genes between groups
    sc.tl.rank_genes_groups(adata, mask_var=avail_mgs, groupby='louvain')
    
    
    ### Construct a dataframe with just genes relevant for assigning cell type
    # After rank_genes_groups, each row in `names` corresponds to each gene, sorted by rank (columns are groups),
    # and each column in `logfoldchanges` correspond to the LFCs of each gene from `names.`
    # This is painful to work with, so instead make a rectangular matrix with columns being groups and rows being
    # only the relevant genes needed to assign cell types.
    
    # n_genes x n_groups dataframe
    gene_group_df = pd.DataFrame(adata.uns["rank_genes_groups"]["names"])
    
    # Retrieve indices for marker genes for each cluster
    # `gene_col` is also the louvain cluster
    gene_row, gene_col = np.where(gene_group_df.isin(avail_mgs))
    
    name_mat = adata.uns["rank_genes_groups"]["names"]
    lfc_mat = adata.uns["rank_genes_groups"]["logfoldchanges"]
    
    mg_group_dict = {
        "cluster": [],
        "gene": [],
        "lfc": []
    }
    
    for i, j in zip(gene_row, gene_col):
        mg_group_dict["cluster"].append(j)
        mg_group_dict["gene"].append(name_mat[i][j])
        mg_group_dict["lfc"].append(lfc_mat[i][j])
    
    mg_group_df = pd.DataFrame.from_dict(mg_group_dict)
    mg_group_df = mg_group_df.pivot(index="cluster", columns="gene").T
    
    # Reformat df after pivot
    mg_group_df = mg_group_df.reset_index()
    mg_group_df.index = mg_group_df["gene"]
    mg_group_df = mg_group_df.drop(["level_0", "gene"], axis=1)
    
    
    cts = list(ct_markers.keys())
    ct_marker_genes = list(ct_markers.values())
    
    # Loop through each group...
    assigned_ct = []
    for i in range(mg_group_df.shape[1]):
        group_lfcs = mg_group_df.iloc[:,i]
    
        # And loop through each cell-type...
        sum_lfcs = []
        for j in range(len(cts)):
            # To sum the log-fold changes over all marker genes in a set
            sum_lfc = group_lfcs.loc[ct_marker_genes[j]].sum()
            sum_lfcs.append(sum_lfc)
    
        # Whichever group of genes associated with a cell-type is highest is what cell type is assigned.
        max_sum = np.argmax(sum_lfcs)
    
        # Enables assignment of "misc" cells
        print(max(sum_lfcs))
        if max(sum_lfcs) < 0:
            assigned_ct.append("Misc.")
        else:
            assigned_ct.append(cts[max_sum])
    
    
    # Map celltype onto louvain groups
    group_ct_map = pd.DataFrame.from_dict({
        "group": np.arange(len(assigned_ct)),
        "celltype": assigned_ct
    })
    
    louvain_groups = pd.DataFrame(adata.obs["louvain"].values.astype(int))
    louvain_groups.columns = ["group"]
    louvain_groups = louvain_groups.merge(group_ct_map)
    
    
    ### Plotting UMAP clusters: Louvain and assigned cell-types
    fig, axes = plt.subplots(1, 2, figsize=(15, 6), layout="constrained")
    
    # Louvain clusters
    axes[0].scatter(umap_x[:,0], umap_x[:,1], c=list(adata.obs["louvain"].values.astype(int)), s=2, alpha=0.8)
    
    u_cts = louvain_groups["celltype"].unique()
    for i in range(len(u_cts)):
        subset_x = np.where(louvain_groups["celltype"] == u_cts[i])[0]
        axes[1].scatter(umap_x[subset_x,0], umap_x[subset_x,1], s=2, alpha=0.8, label=u_cts[i])
    
    fig.supxlabel("UMAP 1")
    fig.supylabel("UMAP 2")
    fig.suptitle(srr_id)
    
    # axes[1].legend(loc="lower center", ncols=3, bbox_to_anchor=(0.5, -0.4))
    axes[1].legend(loc="center right", bbox_to_anchor=(1.45, 0.5))
    
    plt.grid(False)
    plt.savefig(f"../Figures/{srr_id}.jpg", dpi=600, bbox_inches="tight")
    
    
    ### Plotting RNA velocity vectors
    fig, ax = plt.subplots(figsize=(6, 4), dpi=200)
    
    for i in range(len(u_cts)):
        ct_i = np.where(louvain_groups["celltype"] == u_cts[i])[0]
        ax.scatter(umap_x[ct_i,0], umap_x[ct_i,1], s=1, label=u_cts[i], alpha=0.8)
    
    scv.tl.velocity_graph(adata)
    scv.pl.velocity_embedding_stream(adata, basis='umap', ax=ax, smooth=0.25, color="white", show=False)
    
    ax.legend(loc="lower center", ncols=3, bbox_to_anchor=(0.5, -0.2))
    
    plt.savefig(f"../Figures/{srr_id}_velo.jpg", dpi=300, bbox_inches="tight")
    plt.show()
    
    
    ### Save anndata object
    adata.write(f"/nobackup1/ericjlu/{srr_id}.h5ad")
