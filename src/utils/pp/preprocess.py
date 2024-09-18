import scanpy as sc
import numpy as np
from typing import List
import sys
import os
import pandas as pd
import warnings
from tqdm.autonotebook import tqdm
import logging


# Warnings to ignore throughout
sc.settings.verbosity = 0
warnings.simplefilter("ignore", UserWarning)
logging.raiseExceptions = False


# Class to suppress output from inside functions (e.g., PhenoGraph)
class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, "w")

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout


def preprocess(
    adata: sc.AnnData,
    show_progress: bool = True,
    **kwargs,
):
    # Set all arguments below with defaults
    defaults = {
        'hvgs': None,
        'hvgs__n': 2000,
        'hvgs__whitelist': None,
        'pca__var_explained': 0.67,
        'pca__max_comps': 1000,
        'cluster': True,
        'cluster__k': 30,
        'neighbors': True,
        'neighbors__k': 30,
        'umap': True,
        'umap__random_state': 1,
        'impute': True,
        'impute__k': 5,
        'impute__t': 3,
    }
    pp_args = defaults
    pp_args.update(kwargs)

    n_steps = 4
    for key in ['cluster', 'neighbors', 'umap', 'impute']:
        if bool(pp_args[key]): n_steps += 1
    if show_progress: pbar = tqdm(total=n_steps)

    # Median library-size normalization
    if show_progress: pbar.set_description('Normalizing')
    adata.layers['median'] = adata.layers['raw'].copy()
    sc.pp.normalize_total(adata, layer='median')
    if show_progress: pbar.update(1)

    # Log-transformation (natural log, pseudocount of 1)
    if show_progress: pbar.set_description('Log-transforming')
    adata.layers['log'] = adata.layers['median'].copy()
    sc.pp.log1p(adata, layer='log')
    if show_progress: pbar.update(1)

    # Set HVGs
    if show_progress: pbar.set_description('Setting HVGs')
    hvgs = pp_args['hvgs']
    if hvgs is None:
        hvgs = get_hvgs(adata, pp_args['hvgs__n'], pp_args['hvgs__whitelist'])
    adata.var['highly_variable'] = adata.var.index.isin(hvgs)
    n_hvgs = adata.var['highly_variable'].sum()
    if show_progress: pbar.update(1)

    # PCA
    if show_progress: pbar.set_description('Running PCA')
    adata.X = adata.layers['log']
    n_comps = min(pp_args['pca__max_comps'], *adata.shape, n_hvgs) - 2
    sc.tl.pca(adata, n_comps=n_comps, use_highly_variable=True)
    X_pca_full = adata.obsm['X_pca'].copy()
    cum_vars = adata.uns['pca']['variance_ratio'].cumsum()
    n_comps = np.argmin(abs(cum_vars - pp_args['pca__var_explained']))
    adata.obsm['X_pca'] = X_pca_full[:, :n_comps]
    if show_progress: pbar.update(1)

    # Cluster with PhenoGraph
    if pp_args['cluster']:
        if show_progress: pbar.set_description('Clustering with PhenoGraph')
        with HiddenPrints():
            communities, _, _ = sc.external.tl.phenograph(
                pd.DataFrame(adata.obsm['X_pca']),
                k=pp_args['cluster__k'],
                nn_method='brute',
                njobs=-1,
            )
        adata.obs['PhenoGraph_clusters'] = pd.Categorical(communities)
        if show_progress: pbar.update(1)

    # Nearest neighbors in PC space
    if pp_args['neighbors']:
        if show_progress: pbar.set_description('Finding nearest neighbors')
        sc.pp.neighbors(
            adata,
            use_rep='X_pca',
            n_neighbors=pp_args['neighbors__k']
        )
        if show_progress: pbar.update(1)

    # UMAP
    if pp_args['umap']:
        if show_progress: pbar.set_description('Calculating UMAP')
        # Default to PAGA with clusters if clustering, else spectral
        init_pos = 'spectral'
        if pp_args['cluster']:
            sc.tl.paga(adata, groups='PhenoGraph_clusters')
            sc.pl.paga(adata, plot=False)
            init_pos="paga"
        sc.tl.umap(adata, random_state=1, init_pos=init_pos)
        if show_progress: pbar.update(1)

    # Impute expression with MAGIC
    if pp_args['impute']:
        if show_progress: pbar.set_description('Imputing expression with MAGIC')
        adata.X = adata.layers["log"]
        with HiddenPrints():
            try:
                adata_magic = sc.external.pp.magic(
                    adata,
                    copy=True,
                    n_pca=n_comps,
                    knn=pp_args['impute__k'],
                    t=pp_args['impute__t'],
                    verbose=False,
                )
            except ValueError as e:
                pass
        adata.layers['imputed'] = adata_magic.X
        if show_progress: pbar.update(1)

    return adata


def get_hvgs(
    adata: sc.AnnData,
    n_hvgs: int,
    whitelist: List[str],
):
    hvgs = sc.pp.highly_variable_genes(
        adata,
        layer='raw',
        n_top_genes=n_hvgs,
        n_bins=1000,
        flavor='seurat_v3',
        inplace=False,
    )
    hvgs['rank'] = hvgs['highly_variable_rank']

    # Remove genes in blacklist from HVGs
    cwd = os.path.dirname(os.path.realpath(__file__))
    ribosomal = pd.read_csv(
        f'{cwd}/assets/ribosomal_genes.gmt',
        sep='\t',
        index_col=0,
        header=None
    ).loc['RIBOSOMAL_GENES_HUMAN'].iloc[1:].dropna().tolist()
    mitochondrial = hvgs.index[hvgs.index.str.startswith('MT-')].tolist()
    blacklist = set(ribosomal).union(mitochondrial)
    blacklist = hvgs.index.intersection(blacklist)
    for gene in blacklist:
        rank = hvgs.loc[gene, 'rank']
        if not np.isnan(rank):
            hvgs.loc[gene, 'rank'] = np.nan  # remove this gene from rank
            below_rank = hvgs['rank'].gt(rank) & ~hvgs['rank'].isna()
            hvgs.loc[below_rank, 'rank'] -= 1  # move everything else up 1

    # Final list is union of HVGs with whitelist
    hvgs = hvgs.dropna().sort_values('rank').iloc[:n_hvgs]
    hvgs = hvgs.index.union(whitelist).tolist()

    return hvgs