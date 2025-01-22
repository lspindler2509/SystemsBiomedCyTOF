import pandas as pd
import numpy as np
from sklearn import preprocessing
from sklearn.mixture import GaussianMixture as GMM
from scipy.stats import norm
from scipy.spatial.distance import pdist
import scipy.cluster.hierarchy as sch
from collections import Counter
from joblib import Parallel, delayed
import multiprocessing

num_cores = multiprocessing.cpu_count()

def get_label(table):
    return table.index, table.columns

def sort_feature(data):
    pds = np.asarray(data)
    Y = pdist(pds)
    ind = sch.leaves_list(sch.linkage(Y))
    return data[ind, :]


def compute_marker_model(df, table, thres):
    """
    Compute 2-Gaussian mixture model for each marker.
    Outputs mean, variance, and weight for each marker.
    """
    def compute_gmm(mk):
        tmp = df[mk].to_numpy()
        tmp = tmp[tmp > thres]
        if tmp.size == 0:
            return mk, None
        gmm = GMM(n_components=2, n_init=10)
        gmm.fit(tmp[:, np.newaxis])
        index = np.argsort(gmm.means_, axis=0)
        return mk, (
            gmm.means_[index].flatten(),
            gmm.covariances_[index].flatten(),
            gmm.weights_[index].flatten(),
        )

    mk_model = Parallel(n_jobs=num_cores)(delayed(compute_gmm)(mk) for mk in table.columns)
    return {mk: model for mk, model in mk_model if model is not None}


def appr_fun(x, xroot, slope):
    y = np.exp((x - xroot) * slope)
    y = y / (1.0 + y)
    return np.squeeze(y)

def compute_paras(mk_model, mk):
    mus, sigmas, ws = mk_model[mk]
    a = (-0.5 / sigmas[0] + 0.5 / sigmas[1])
    b = mus[0] / sigmas[0] - mus[1] / sigmas[1]
    c = (
        0.5 * (-mus[0] ** 2 / sigmas[0] + mus[1] ** 2 / sigmas[1])
        + np.log(ws[0] / ws[1])
        + 0.5 * np.log(sigmas[1] / sigmas[0])
    )
    xroot = (-b - np.sqrt(b ** 2 - 4.0 * a * c)) / (2.0 * a)
    slope = 0.5 * (xroot - mus[0]) / sigmas[0] - 0.5 * (xroot - mus[1]) / sigmas[1]
    return xroot, slope

def score_fun(x, mk_model, mk):
    xroot, slope = compute_paras(mk_model, mk)
    return appr_fun(x, xroot, slope)

def get_score_mat(X, weights, table, y_cluster, mk_model):
    if y_cluster is not None and len(y_cluster) > 0:
        clusters = np.unique(y_cluster)
        if len(clusters) > 0:
            mu = np.vstack(
                [np.mean(X[y_cluster == cl, :], axis=0) for cl in clusters if np.any(y_cluster == cl)]
            )
            score = get_score_mat_main(mu, weights, table, mk_model)
        else:
            raise ValueError("No valid clusters found in y_cluster.")
    else:
        score = get_score_mat_main(X, weights, table, mk_model)
    return score

def get_score_mat_main(X, weights, table, mk_model):
    score = np.zeros((X.shape[0], len(table.index)))
    
    def process_row(i):
        score_tmp = np.ones(X.shape[0])
        for j, mk in enumerate(table.columns):
            if table.loc[i, mk] > 0:
                score_tmp = np.minimum(score_tmp, score_fun(X[:, j], mk_model, mk))
            elif table.loc[i, mk] < 0:
                score_tmp = np.minimum(score_tmp, 1.0 - score_fun(X[:, j], mk_model, mk))
        return score_tmp

    results = Parallel(n_jobs=num_cores)(delayed(process_row)(ct) for ct in table.index)
    for i, result in enumerate(results):
        score[:, i] = result
    return score


def get_unique_index(X, score, table, thres):
    """
    Get unique indices for clusters based on scores.
    """
    ct_score = np.abs(table.values).sum(axis=1)  # Replaced as_matrix() with values
    ct_index = np.zeros((X.shape[0], len(table.index) + 1))

    for i, ct in enumerate(table.index):
        ct_index[:, i] = (score[:, i] > thres) * 1

    for i in np.where(ct_index.sum(axis=1) > 1)[0]:
        mk_tmp = np.where(ct_index[i, :] == 1)[0]
        score_tmp = ct_score[mk_tmp]
        ct_index[i, :] = 0
        if np.sum(score_tmp == score_tmp.max()) == 1:
            ct_index[i, mk_tmp[np.argmax(score_tmp)]] = 1

    ct_index[:, -1] = (score[:, -1] > thres) * 1
    return ct_index


def get_landmarks(X, score, ct_index, idx2ct, obj, thres):
    """
    Extract landmark points for visualization or analysis.
    """
    res_c = {}
    for idx, ct in enumerate(idx2ct):
        if np.any(ct_index[:, idx] == 1):
            item = X[ct_index[:, idx] == 1, :]
        else:
            item = X[score[:, idx] > thres, :]

        if item.shape[0] > 60:
            res = obj.cluster(
                item,
                k=30,
                directed=False,
                prune=False,
                min_cluster_size=10,
                jaccard=True,
                primary_metric="euclidean",
                n_jobs=-1,
                q_tol=1e-3,
            )
            res_c[ct] = return_center(item, res[0])
        elif item.shape[0] > 0:
            res_c[ct] = item.mean(axis=0)[np.newaxis, :]

    return res_c


def return_center(X, labels):
    """
    Compute cluster centers.
    """
    clusters = np.unique(labels)
    centers = []
    for i in clusters:
        centers.append(np.mean(X[labels == i, :], axis=0))
    return np.vstack(centers)


def select_centers(res_center, X, table, mk_model):
    """
    Select representative centers for each cluster.
    """
    res_score = {}
    res_select = {}

    for ct, item in res_center.items():
        if item.shape[0] > 1:
            score = np.ones(item.shape[0])
            for j, mk in enumerate(table.columns):
                if table.loc[ct, mk] == 1:
                    score = np.minimum(score, score_fun(item[:, j], mk_model, mk))
                elif table.loc[ct, mk] == -1:
                    score = np.minimum(score, 1.0 - score_fun(item[:, j], mk_model, mk))

            res_score[ct] = score
            res_select[ct] = res_center[ct][np.argmax(score), :][np.newaxis, :]
        else:
            res_select[ct] = res_center[ct].copy()
    return res_score, res_select


def output_feature_matrix(res, label=None):
    """
    Convert a dictionary structure into a feature matrix.
    """
    ct_mat = []
    feature_mat = []
    if label is None:
        label = res.keys()

    for key in label:
        if key in res:
            item = res[key]
            if item.shape[0] > 0:
                feature_mat.append(item)
                ct_mat += [key] * item.shape[0]

    feature_mat = np.vstack(feature_mat)
    return feature_mat, ct_mat




# Compute scores using marker-specific thresholds
def get_score_mat_with_thresholds(X, weights, table, y_cluster, mk_model, thresholds):
    """
    Computes the score matrix while using marker-specific thresholds.
    """
    if y_cluster is not None and len(y_cluster) > 0:
        clusters = np.unique(y_cluster)
        if len(clusters) > 0:
            mu = np.vstack(
                [np.mean(X[y_cluster == cl, :], axis=0) for cl in clusters if np.any(y_cluster == cl)]
            )
            score = get_score_mat_main_with_thresholds(mu, weights, table, mk_model, thresholds)
        else:
            raise ValueError("No valid clusters found in y_cluster.")
    else:
        score = get_score_mat_main_with_thresholds(X, weights, table, mk_model, thresholds)
    return score

def get_score_mat_main_with_thresholds(X, weights, table, mk_model, thresholds):
    score = np.zeros((X.shape[0], len(table.index)))
    
    def process_row(i):
        score_tmp = np.ones(X.shape[0])
        for j, mk in enumerate(table.columns):
            if table.loc[i, mk] > 0:
                score_tmp = np.minimum(score_tmp, appr_fun(X[:, j], *compute_paras_with_threshold(mk_model, mk, thresholds[mk])))
            elif table.loc[i, mk] < 0:
                score_tmp = np.minimum(score_tmp, 1.0 - appr_fun(X[:, j], *compute_paras_with_threshold(mk_model, mk, thresholds[mk])))
        return score_tmp

    results = Parallel(n_jobs=num_cores)(delayed(process_row)(ct) for ct in table.index)
    for i, result in enumerate(results):
        score[:, i] = result
    return score

def compute_paras_with_threshold(mk_model, mk, threshold):
    """
    Compute parameters for score function considering marker-specific thresholds.
    """
    mus, sigmas, ws = mk_model[mk]
    a = (-0.5 / sigmas[0] + 0.5 / sigmas[1])
    b = mus[0] / sigmas[0] - mus[1] / sigmas[1]
    c = (
        0.5 * (-mus[0] ** 2 / sigmas[0] + mus[1] ** 2 / sigmas[1])
        + np.log(ws[0] / ws[1])
        + 0.5 * np.log(sigmas[1] / sigmas[0])
    )
    xroot = (-b - np.sqrt(b ** 2 - 4.0 * a * (c - threshold))) / (2.0 * a)
    slope = 0.5 * (xroot - mus[0]) / sigmas[0] - 0.5 * (xroot - mus[1]) / sigmas[1]
    return xroot, slope
