import pandas as pd
import numpy as np
from random import shuffle
import scipy.stats as stats

def cal_flux(data, elements, significant, overall_mean):
    enrichment = []
    num_cluster = len(np.unique(data["assignment"].tolist()))
    for c in np.unique(data["assignment"].tolist()):
        data_run = data.copy()
        cluster = data_run[data_run["assignment"] == c]
        data_run.loc[
            data_run[data_run["assignment"] != c].index, "assignment"
        ] = "R"
        for el in elements:
            if cluster[el].mean() > overall_mean[el]:
                diff = cluster[el].mean() - overall_mean[el]
                contingency_table = pd.crosstab(
                    data_run[el], data_run["assignment"]
                )
                p_val = stats.fisher_exact(contingency_table)[1]
                if significant:
                    if p_val < (0.05 / (num_cluster * 6)):
                        enrichment.append([c, el, np.round(diff, 3)])
                else:
                    enrichment.append([c, el, np.round(diff, 3)])
    sign_diff = pd.DataFrame(np.array(enrichment))
    sign_diff.columns = ["cluster", "element", "diff"]
    if '.' in sign_diff["element"].iloc[0]:
        sign_diff["element"] = sign_diff["element"].str.split(".", expand=True)[1]
    elif '_' in sign_diff["element"].iloc[0]: 
        sign_diff["element"] = sign_diff["element"].str.split("_", expand=True)[1]
    sign_diff["cluster"]  = sign_diff.cluster.astype(float).astype(int)
    names = sign_diff.cluster.unique().tolist() + sign_diff.element.unique().tolist()
    sign_diff_expand = pd.DataFrame(columns=names, index=names)
    for c in sign_diff.cluster.unique():
        for el in sign_diff.element.unique():
            diff = (
                sign_diff[(sign_diff["cluster"] == c) & (sign_diff["element"] == el)][
                    "diff"
                ]
                .astype(float)
                .values
            )
            if len(diff) != 0:
                sign_diff_expand.loc[c, el] = diff[0]

    sign_diff_expand = sign_diff_expand.fillna(0)
    sign_diff_expand = sign_diff_expand.div(
        sign_diff_expand.sum(axis=1), axis=0
    ).fillna(0)
    return sign_diff_expand, names