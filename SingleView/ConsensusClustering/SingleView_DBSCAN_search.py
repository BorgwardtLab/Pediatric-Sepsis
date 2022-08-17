import pandas as pd
pd.set_option("display.precision", 3)
import os
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from ConsensusClusteringSingleView import ConsensusCluster
import scipy.stats as sps
import copy
from tqdm import tqdm
from sklearn.manifold import TSNE
from sklearn.metrics import silhouette_score
from sklearn.cluster import DBSCAN
from tqdm import tqdm
import argparse
from pathlib import Path
from itertools import product
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--view", type=str, default="clinical")
    args = parser.parse_args()
    data_path = "data/"
    score_path = "{}/Clustering_silhouette/".format(data_path)
    cdf_path = "{}/CDF plots/".format(data_path)
    tsne_path = "{}/TSNEplots/".format(data_path)
    KCC_path = "{}/KCC/".format(data_path)
    f_stat_path = '{}/F_stat/'.format(data_path)
        
    KCC_spaces = np.arange(3, 11)
    n_states = np.arange(0, 1)
    epss = np.linspace(0.1, 10, 50)
    minss = np.arange(2, 22)
    DBSCAN_grids = list(product(n_states, epss, minss))
    perplexitys = [30]
    random_states = np.arange(0, 5)
    TSNE_grids = list(product(perplexitys, random_states))
    perf = pd.DataFrame(
        columns=[
            "View",
            "KCC_space",
            "TSNE",
            "perplexity",
            "random_state",
            "eps",
            "min_samples",
            "n_state",
            "silhouette_score",
            "DBSCAN_clusters",
            "outlier_ratio",
        ]
    )
    j = 0
    for KCC_space in range(3, 11):
        data = pd.read_csv(
            "{}/KCC_Cov_AgeSexEth_{}_NormalImputation_K{}.csv".format(
                KCC_path, args.view, KCC_space
            ),
            header=None,
        )
        assert data.shape == (387, KCC_space)
        for DBSCAN_param in tqdm(DBSCAN_grids):
            n_state, eps, mins = DBSCAN_param
            np.random.seed(n_state)
            cluster = DBSCAN(eps=eps, min_samples=mins)
            for TSNE_param in TSNE_grids:
                perplexity, random_state = TSNE_param
                tsne = TSNE(
                        perplexity=perplexity, random_state=random_state, n_components=2
                    )
                X_emb = tsne.fit_transform(data.values)
                cluster.fit(X_emb)
                outlier_ratio = (cluster.labels_ == -1).sum() / len(data)
                if len(np.unique(cluster.labels_[cluster.labels_ != -1])) > 1:
                    score = silhouette_score(
                        X_emb[cluster.labels_ != -1],
                        cluster.labels_[cluster.labels_ != -1],
                    )
                    DBSCAN_clusters = len(
                        np.unique(cluster.labels_[cluster.labels_ != -1])
                    )
                    entry = [
                        args.view,
                        KCC_space,
                        True,
                        perplexity,
                        random_state,
                        eps,
                        mins,
                        n_state,
                        score,
                        DBSCAN_clusters,
                        outlier_ratio,
                    ]
                    perf.loc[j] = entry
                    j += 1
        perf.to_csv("{}/DBSCAN_{}_view.csv".format(score_path, args.view))
if __name__ == "__main__":
    main()