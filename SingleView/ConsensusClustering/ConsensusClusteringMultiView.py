#####################################################################
# Modified from https://github.com/ZigaSajovic/Consensus_Clustering #
#####################################################################

import numpy as np
from itertools import combinations
import bisect
from random import sample


class ConsensusCluster:
    """
    Implementation of Consensus clustering, following the paper
    https://link.springer.com/content/pdf/10.1023%2FA%3A1023949509487.pdf
    Args:
      * cluster -> clustering class
      * NOTE: the class is to be instantiated with parameter `n_clusters`,
        and possess a `fit_predict` method, which is invoked on data.
      * L -> smallest number of clusters to try
      * K -> biggest number of clusters to try
      * H -> number of resamplings for each cluster number
      * resample_proportion -> percentage to sample
      * Mk -> consensus matrices for each k (shape =(K,data.shape[0],data.shape[0]))
              (NOTE: every consensus matrix is retained, like specified in the paper)
      * Ak -> area under CDF for each number of clusters
              (see paper: section 3.3.1. Consensus distribution.)
      * deltaK -> changes in areas under CDF
              (see paper: section 3.3.1. Consensus distribution.)
      * self.bestK -> number of clusters that was found to be best
    """

    def __init__(self, cluster, L, K, H, resample_proportion=0.5):
        assert 0 <= resample_proportion <= 1, "proportion has to be between 0 and 1"
        self.cluster_ = cluster
        self.resample_proportion_ = resample_proportion
        self.L_ = L
        self.K_ = K
        self.H_ = H
        self.Mk = None
        self.Ak = None
        self.deltaK = None
        self.bestK = None

    def _internal_resample(self, Xs, proportion):
        """
        Args:
          * data -> (examples,attributes) format
          * proportion -> percentage to sample
        """
        resampled_indices = sample(
            range(Xs[0].shape[0]), int(Xs[0].shape[0] * proportion)
        )
        return (
            resampled_indices,
            Xs[0].values[resampled_indices, :],
            Xs[1].values[resampled_indices, :],
        )

    def fit(self, Xs, verbose=False):
        """
        Fits a consensus matrix for each number of clusters
        Args:
          * data -> (examples,attributes) format
          * verbose -> should print or not
        """
        Mk = np.zeros((self.K_ - self.L_, Xs[0].shape[0], Xs[0].shape[0]))
        Is = np.zeros((Xs[0].shape[0],) * 2)
        for k in range(self.L_, self.K_):  # for each number of clusters
            i_ = k - self.L_
            if verbose:
                print("At k = %d, aka. iteration = %d" % (k, i_))
            for h in range(self.H_):  # resample H times
                if verbose:
                    print("\tAt resampling h = %d, (k = %d)" % (h, k))
                (
                    resampled_indices,
                    resample_Clinical,
                    resample_NPX,
                ) = self._internal_resample(Xs, self.resample_proportion_)
                Xs_sampled = [resample_Clinical, resample_NPX]
                Mh = self.cluster_(
                    n_clusters=k, random_state=10, n_init=10
                ).fit_predict(Xs_sampled)
                # find indexes of elements from same clusters with bisection
                # on sorted array => this is more efficient than brute force search
                index_mapping = np.array((Mh, resampled_indices)).T
                # index_mapping = index_mapping[index_mapping[:, 0].argsort()]
                clusts_ = index_mapping[:, 0]
                id_clusts = index_mapping[:, 1]
                for i in range(k):
                    is_ = id_clusts[clusts_ == i]
                    ids_ = np.array(list(combinations(is_, 2))).T
                    assert np.isin(is_, resampled_indices).all()
                    if ids_.size != 0:
                        Mk[i_, ids_[0], ids_[1]] += 1
                        Mk[i_, ids_[1], ids_[0]] += 1
                ids_2 = np.array(list(combinations(resampled_indices, 2))).T
                Is[ids_2[0], ids_2[1]] += 1
                Is[ids_2[1], ids_2[0]] += 1
            assert (Mk[i_] == Mk[i_].T).all()
            assert (Is[i_] == Is[i_].T).all()
            Mk[i_] /= Is + 1e-8
            Mk[i_, range(Xs[0].shape[0]), range(Xs[0].shape[0])] = 1  # always with self
            Is.fill(0)  # reset counter
            assert (Mk[i_] == Mk[i_].T).all()
            assert np.max(Mk[i_]) <= 1
        self.Mk = Mk
        # fits areas under the CDFs
        self.Ak = np.zeros(self.K_ - self.L_)
        for i, m in enumerate(Mk):
            hist, bins = np.histogram(m.ravel(), density=True)
            self.Ak[i] = np.sum(
                h * (b - a) for b, a, h in zip(bins[1:], bins[:-1], np.cumsum(hist))
            )
        # fits differences between areas under CDFs
        self.deltaK = np.array(
            [
                (Ab - Aa) / Aa if i > 2 else Aa
                for Ab, Aa, i in zip(
                    self.Ak[1:], self.Ak[:-1], range(self.L_, self.K_ - 1)
                )
            ]
        )
        self.bestK = (
            np.argmax(self.deltaK) + self.L_ if self.deltaK.size > 0 else self.L_
        )

    def predict(self):
        """
        Predicts on the consensus matrix, for best found cluster number
        """
        assert self.Mk is not None, "First run fit"
        return self.cluster_(n_clusters=self.bestK).fit_predict(
            1 - self.Mk[self.bestK - self.L_]
        )

    def predict_data(self, Xs):
        """
        Predicts on the data, for best found cluster number
        Args:
          * data -> (examples,attributes) format
        """
        assert self.Mk is not None, "First run fit"
        return self.cluster_(n_clusters=self.bestK).fit_predict(Xs)
