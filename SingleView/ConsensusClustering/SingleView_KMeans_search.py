import pandas as pd

pd.set_option("display.precision", 3)
import os
import warnings

warnings.filterwarnings("ignore")
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
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
from ConsensusClusteringSingleView import ConsensusCluster

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

    perplexitys = [5, 10, 20, 30, 50]
    random_states = np.arange(0, 10)
    TSNE_grids = list(product(perplexitys, random_states))

    perf = pd.DataFrame(
        columns=[
            "View",
            "KCC_space",
            "TSNE",
            "perplexity",
            "random_state",
            "silhouette_score"
        ]
    )
    j = 1
    for KCC_space in tqdm(range(3, 11)):
        data = pd.read_csv(
            "{}/KCC_Cov_AgeSexEth_{}_NormalImputation_K{}.csv".format(
                KCC_path, args.view, KCC_space
            ),
            header=None,
        )
        assert data.shape == (387, KCC_space)
        cons = ConsensusCluster(KMeans, KCC_space, KCC_space + 1, 100, resample_proportion=0.8)
        for use_TSNE in [True, False]:
            if use_TSNE:
                for TSNE_param in TSNE_grids:
                    perplexity, random_state = TSNE_param
                    tsne = TSNE(
                        perplexity=perplexity, random_state=random_state, n_components=2
                    )
                    X_emb = tsne.fit_transform(data.values)
                    assert X_emb.shape == (387, 2)
                    cons.fit(pd.DataFrame(X_emb))
                    assignment = cons.predict_data(X_emb) + 1
                    score = silhouette_score(X_emb, assignment)
                    entry = [args.view, KCC_space, use_TSNE, perplexity, random_state, score]
                    perf.loc[j] = entry
                    j += 1
            else:
                cons.fit(data)
                assignment = cons.predict_data(data.values) + 1
                score = silhouette_score(data.values, assignment)    
                entry = [args.view, KCC_space, use_TSNE, "NA", "NA", score]
                perf.loc[j] = entry
                j += 1

        perf.to_csv("{}/ConsensusKMeans_{}_view.csv".format(score_path, args.view))

if __name__ == "__main__":
    main()
        
       

