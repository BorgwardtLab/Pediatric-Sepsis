import pandas as pd

pd.set_option("display.precision", 3)
import os
import warnings

warnings.filterwarnings("ignore")
from sklearn.cluster import SpectralClustering
from mvlearn.cluster import MultiviewKMeans
from sklearn.cluster import KMeans
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from IPython.display import display
import sys
import math
import seaborn as sns
import argparse

sns.set_style("white")
from ConsensusClusteringSingleView import ConsensusCluster
import scipy.stats as sps
import copy
from tqdm import tqdm
from sklearn.manifold import TSNE
import matplotlib.cm as cm
import matplotlib.lines as mlines
from sklearn.metrics import silhouette_score


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--view", type=str, default="clinical")
    args = parser.parse_args()
    PCA = "_PCA_standardized"
    data_path = "data{}/".format(PCA)
    score_path = "Clustering_silhouette{}/".format(PCA)
    cdf_path = "CDF plots{}/".format(PCA)
    tsne_path = "TSNEplots{}/".format(PCA)

    perf = pd.DataFrame(columns=["view", "KCC_space", "KMeans_cluster", "silhouette"])

    view = args.view
    i = 0
    for kmeans_cluster in tqdm(range(2, 21)):
        Mks = []
        Aks = []
        assignments = pd.DataFrame(columns=range(2, 21))
        for kcc_space in range(2, 21):
            data = pd.read_csv(
                "{}/KCC_Cov_AgeSexEth_{}_PCA_NormalImputation_PCA95%_K{}.csv".format(
                    data_path, view, kcc_space
                ),
                header=None,
            )

            assert data.shape == (387, kcc_space)
            cons = ConsensusCluster(
                KMeans, kmeans_cluster, kmeans_cluster + 1, 100, resample_proportion=0.8
            )
            cons.fit(data)
            assignments[kcc_space] = cons.predict_data(data) + 1
            assert cons.Mk.shape[0] == 1
            Mks.append(cons.Mk[0])
            score = silhouette_score(data.values, assignments[kcc_space].tolist())
            entry = [view, kcc_space, kmeans_cluster, score]
            perf.loc[i] = entry
            i += 1

        assignments.to_csv(
            "{}/ConsensusAssignments_Cov_AgeSexEth_{}_KMeans_{}_NormalImputation_PCA95%.csv".format(
                data_path, view, kmeans_cluster
            )
        )

        Aks = []
        plt.figure(figsize=(20, 15))
        for kcc_space in range(2, 21):
            hist, bins = np.histogram(Mks[kcc_space - 2].ravel(), bins=50)
            pdf = hist / sum(hist)
            cdf = np.cumsum(pdf)
            acdf = np.sum(h * (b - a) for b, a, h in zip(bins[1:], bins[:-1], cdf))
            Aks.append(acdf)
            plt.plot(bins[1:], cdf, label=kcc_space)

        plt.legend(fontsize=16)
        plt.title("CDF: {} view".format(view), fontsize=25)
        plt.ylim(-0.05, 1.05)
        plt.tight_layout()
        plt.savefig(
            "{}/Cov_AgeSexEth_{}_KMeans_{}_NormalImputation_PCA95%.png".format(
                cdf_path, view, kmeans_cluster
            ),
            dpi=300,
        )
        plt.close()
    perf.to_csv("{}/KCC_ConsensusKMeans_{}_view.csv".format(score_path, view))


if __name__ == "__main__":
    main()
