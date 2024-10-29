import pandas as pd
import scipy as sp
from typing import List, Tuple


def get_binary_association(
    enrichments: pd.DataFrame,
    annotations: pd.DataFrame,
    module: str,
    column: str,
    negative_vals: List,
):
    idx = enrichments[module].dropna().index
    idx = idx.intersection(annotations[column].dropna().index)
    mask = annotations.loc[idx, column].isin(negative_vals)
    r, p = sp.stats.ranksums(
        enrichments.loc[idx[~mask], module],
        enrichments.loc[idx[mask], module],
    )
    return r, p


def get_continuous_association(
    enrichments: pd.DataFrame,
    annotations: pd.DataFrame,
    module: str,
    column: str,
):
    idx = enrichments[module].dropna().index
    idx = idx.intersection(annotations[column].dropna().index)
    r, p = sp.stats.pearsonr(enrichments[module])
    return r, p


def get_associations(
    enrichments: pd.DataFrame,
    annotations: pd.DataFrame,
    modules: List,
    binary_columns: dict,
    continuous_columns: List[str],
):
    associations = pd.DataFrame(
        index=modules,
        columns=list(binary_columns.keys()) + continuous_columns,
        dtype=float,
    )
    pvals = associations.copy()

    for module in modules:
        for column in continuous_columns:
            r, p = get_continuous_association(
                enrichments, annotations, module, column, 
            )
            associations.loc[module, column] = r
            pvals.loc[module, column] = p
        
        for column, negative in binary_columns.items():
            if type(negative) != list: negative = [negative]
            r, p = get_binary_association(
                enrichments, annotations, module, column, negative,
            )
            associations.loc[module, column] = r
            pvals.loc[module, column] = p

    return associations, pvals
