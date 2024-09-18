import scipy as sp
import numpy as np
import scanpy as sc
import pandas as pd
from sklearn.neighbors import NearestNeighbors
from phenograph.classify import random_walk_probabilities
from sklearn.linear_model import LinearRegression
import harmony.core


# Get affinity matrix for subset of AnnData (not augmented)
def get_affinity_matrix(
    adata: sc.AnnData,
    n_neighbors: int, 
    metric: str,
):
    # Get nearest neighbors
    sc.pp.neighbors(
        adata,
        use_rep='X_pca',
        n_neighbors=n_neighbors,
        metric=metric,
    )
    dists = adata.obsp['distances']

    # Get the (adaptive) kth nearest neighbor for each cell
    adaptive_k = int(np.floor(n_neighbors/3))
    adaptive_std = np.zeros(adata.shape[0])
    for i, chunk in enumerate(np.split(dists.data, dists.indptr[1:-1])):
        adaptive_std[i] = np.sort(chunk)[adaptive_k-1]

    # Normalize dists and return diffusion kernel
    i, j, d = sp.sparse.find(dists)
    d /= adaptive_std[i] # normalize dists by the (adaptive) kth nearest neighbor
    W = sp.sparse.csr_matrix(
        (np.exp(-d), (i,j)),
        shape=adata.obsp['distances'].shape,
    )
    kernel = W + W.T

    return kernel, pd.Series(adaptive_std, index=adata.obs.index)


# Function to construct mutually nearest neighbors bewteen two datasets
def construct_mnn(
    t1_cells,
    t2_cells,
    data_df: pd.DataFrame,
    n_neighbors: int,
    metric: str,
    n_jobs=-2,
):
    nbrs = NearestNeighbors(
        n_neighbors=n_neighbors,
        metric=metric,
        n_jobs=n_jobs
    )
    t1_data = data_df.loc[t1_cells, :].values
    t2_data = data_df.loc[t2_cells, :].values
    # Dataset 1 neighbors
    nbrs.fit(t1_data)
    t1_nbrs = nbrs.kneighbors_graph(t2_data, mode='distance')

    # Dataset 2 neighbors
    nbrs.fit(t2_data)
    t2_nbrs = nbrs.kneighbors_graph(t1_data, mode='distance')

    # Mututally nearest neighbors
    mnn = t2_nbrs.multiply(t1_nbrs.T)
    mnn = mnn.sqrt()
    return mnn


# From Harmony
def mnn_ka_distances(mnn, n_neighbors):
    # Function to find distance kth neighbor in the mutual nearest neighbor matrix
    ka = int(n_neighbors / 3)
    ka_dists = np.repeat(None, mnn.shape[0])
    x, y, z = sp.sparse.find(mnn)
    rows=pd.Series(x).value_counts()
    for r in rows.index[rows >= ka]:
        ka_dists[r] = np.sort(z[x==r])[ka - 1]
    return ka_dists


# From Harmony
def mnn_scaling_factors(mnn_ka_dists, scaling_factors):
    cells = mnn_ka_dists.index[~mnn_ka_dists.isnull()]
    # Linear model fit
    x = scaling_factors[cells]
    y = mnn_ka_dists[cells]
    lm = LinearRegression()
    lm.fit(x.values.reshape(-1, 1), y.values.reshape(-1, 1))
    # Predict
    x = scaling_factors[mnn_ka_dists.index]
    vals = np.ravel(lm.predict(x.values.reshape(-1, 1)))
    mnn_scaling_factors = pd.Series(vals, index=mnn_ka_dists.index)
    return mnn_scaling_factors


def get_mnn_affinity_function(
    index_1,
    index_2,
    scaling_factors,
    adata,
    n_neighbors: int,
    metric: str,
):
    sorted_index = index_1.append(index_2)

    # Construct MNN between groups
    mnn = construct_mnn(
        index_1, index_2, 
        pd.DataFrame(adata.obsm['X_pca'], index=adata.obs.index), 
        n_neighbors,
        metric=metric
    )
    # MNN adaptive distances
    ka_dists = pd.Series(0.0, index=sorted_index) # distance to (adaptive) kth neighbor
    ka_dists[index_1] = mnn_ka_distances(mnn, n_neighbors)
    ka_dists[index_2] = mnn_ka_distances(mnn.T, n_neighbors)

    # MNN scaling factors
    mnn_sf = pd.Series(0.0, index=sorted_index)
    mnn_sf[index_1] = mnn_scaling_factors(
        ka_dists[index_1], scaling_factors,
    )
    mnn_sf[index_2] = mnn_scaling_factors(
        ka_dists[index_2], scaling_factors,
    )
    # MNN affinity matrix
    mnn_aff = harmony.core._mnn_affinity(
        mnn, mnn_sf,
        np.where(adata.obs.index.isin(index_1))[0][0],
        np.where(adata.obs.index.isin(index_2))[0][0],
        device='cpu',
    )
    return mnn_aff


# Calculate augmented affinity matrix between labeled and unlabeled datasets
def get_augmented_affinity_matrix(
    adata: sc.AnnData,
    label_column: str,
    knn: int = 30,
    mnn: int = 60,
):
    # Assume unlabeled samples have NaN in label column
    is_labeled = ~adata.obs[label_column].isna()

    # Affinity matrix for unlabeled data
    aff_unl, sf_unl = get_affinity_matrix(
        adata[~is_labeled],
        n_neighbors=knn,
        metric='euclidean',
    )
    # Affinity matrix for labeled data
    aff_lbl, sf_lbl = get_affinity_matrix(
        adata[is_labeled],
        n_neighbors=knn,
        metric='euclidean',
    )

    # Affinity matrix between labeled and unlabeled data
    sf_cmb = pd.concat([sf_unl, sf_lbl])
    adata = adata[sf_cmb.index].copy() # reorder AnnData
    mnn_aff = get_mnn_affinity_function(
        sf_unl.index,
        sf_lbl.index,
        sf_cmb,
        adata,
        n_neighbors=mnn,
        metric='cosine',
    )

    # Combine affinity matrices
    # Expand shape of affinity matrices to shape of combined AnnData
    shape = [adata.obs.shape[0]]*2 # cells x cells

    # In-vitro affinity matrix
    i, j, v = sp.sparse.find(aff_unl)
    aff_unl_full = sp.sparse.csr_matrix((v, (i, j)), shape=shape)

    # In-vivo affinity matrix
    i, j, v = sp.sparse.find(aff_lbl)
    offset = aff_unl.shape[0] 
    aff_lbl_full = sp.sparse.csr_matrix((v, (i+offset, j+offset)), shape=shape)

    # Combine all affinity matrices into symmetric matrix
    comb_aff = aff_unl_full + aff_lbl_full + mnn_aff + mnn_aff.T

    return comb_aff


def transfer_labels(
    adata,
    label_column: str,
    knn: int = 30,
    mnn: int = 60,
):
    # Assume unlabeled samples have NaN in label column
    is_labeled = ~adata.obs[label_column].isna()

    # Get labels for known and unknown data
    labels = adata.obs.loc[is_labeled, label_column].values
    lbl_codes, lbl_uniques = pd.factorize(labels)
    lbl_codes = np.append(
        np.zeros((np.sum(~is_labeled),), dtype=int),  # unlabeled is zero
        lbl_codes + 1
    )
    
    # Map labels using PhenoGraph classify and affinity matrix
    A = get_augmented_affinity_matrix(adata, label_column, knn, mnn)
    P = random_walk_probabilities(A, lbl_codes)
    c = np.argmax(P, axis=1)
    
    # Annotate AnnData
    c_map = dict(enumerate(lbl_uniques))
    adata.obs.loc[~is_labeled, label_column] = pd.Series(c).map(c_map).values
    for key, val in c_map.items():
        adata.obs[f'P({val})'] = 0.
        adata.obs.loc[is_labeled, f'P({val})'] = 1.
        adata.obs.loc[~is_labeled, f'P({val})'] = P.T[key]
