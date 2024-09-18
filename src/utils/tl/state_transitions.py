import pandas as pd
import scanpy as sc
import itertools
from typing import List
import numpy as np
import scipy as sp
import tqdm


def get_feature_overlaps(
    adata: sc.AnnData,
    feature_columns: List[str],
    qtl_threshold: float = 0.75,
) -> pd.DataFrame:
    # Count co-occurrence of features in cells
    thresholds = adata.obs[feature_columns].quantile(qtl_threshold)
    counts = adata.obs[feature_columns].ge(thresholds)
    overlaps = pd.DataFrame(dtype=float)
    for f1 in feature_columns:
        for f2 in feature_columns:
            overlaps.loc[f1, f2] = counts[[f1, f2]].all(axis=1).sum()
    # Frequencies of co-occurrence of features in cells
    overlaps = overlaps.divide(overlaps.max(0), 0)
    return overlaps


def test_ratios(
    df: pd.DataFrame,
    groupby: List[str],
    n_iters: int = 1000,
):
    # Reference fractions
    fracs = df.groupby(groupby).mean()['High']
    ratios = fracs.xs('Metastasis', level=1) / fracs.xs('Primary', level=1)
    mask = ~ratios.isin([np.inf, 0])
    ratios = ratios.replace(np.inf, ratios[mask].max())
    ratios = ratios.replace(0, ratios[mask].min())

    # Null fractions from random shuffling
    null_ratios = pd.DataFrame(
        index=ratios.index, 
        columns=range(n_iters),
        dtype=float,
    )

    def shuffle(group):
        idx = group.index
        group = group.sample(frac=1)
        group.index = idx
        return group

    for i in tqdm.tqdm(null_ratios.columns):
        null_df = df.copy()
        gb = null_df.groupby('Patient', group_keys=False)
        null_df['High'] = gb['High'].apply(shuffle)
        null_fracs = null_df.groupby(groupby).mean()['High']
        null_ratios[i] = null_fracs.xs('Metastasis', level=1) 
        null_ratios[i] /= null_fracs.xs('Primary', level=1)

    # Remove infinities for stats calculation
    mask = ~null_ratios.isin([np.inf, 0])
    null_ratios = null_ratios.replace(np.inf, null_ratios[mask].max())
    null_ratios = null_ratios.replace(0, null_ratios[mask].min())
    r, p = sp.stats.ranksums(
        np.log10(ratios.values.flatten()),
        np.log10(null_ratios.values.flatten()),
        alternative='greater',
    )

    return r, p
