from cdlib import algorithms
import networkx as nx
import pandas as pd 
import numpy as np
import pickle as pkl 
import plotly.graph_objects as go
from scipy.stats import ttest_ind
from collections import Counter 
import seaborn as sns
import matplotlib.pyplot as plt
import os
from scipy.sparse.csgraph import connected_components
from scipy.sparse import csr_matrix 


def get_cohen_d(x, y):
  """
  Calculates Cohen's d effect size for two independent samples.

  Args:
    x: A list or NumPy array containing the first sample.
    y: A list or NumPy array containing the second sample.

  Returns:
    The calculated Cohen's d value.
  """

  nx = len(x)
  ny = len(y)

  dof = nx + ny - 2

  pooled_std = np.sqrt(((nx-1)*np.std(x, ddof=1)**2 + (ny-1)*np.std(y, ddof=1)**2) / dof)
  std_x = np.std(x, ddof=1)
  std_y = np.std(y, ddof=1)
  return (np.mean(x) - np.mean(y)) / pooled_std, std_x, std_y


def shuffle_symmetric_matrix(matrix):
    idx = matrix.index
    matrix = matrix.to_numpy()
    if not np.allclose(matrix, matrix.T):
        raise ValueError("Input matrix must be symmetric.")

    # Generate a random permutation of indices
    n = matrix.shape[0]
    perm = np.random.permutation(n)

    # Apply the permutation to rows and columns
    shuffled_matrix = matrix[perm][:, perm]
    shuffled_matrix = pd.DataFrame(shuffled_matrix, columns= idx, index=idx)
    return shuffled_matrix


SEQ_PATH = "../pm_clustering/results/proteins/seq_len_dict.pkl"
SEQ_DICT = pkl.load(open(SEQ_PATH, "rb"))
DB_PATH = "../Pocketeome/db_files"
PLOT_DIR = "../pm_clustering/results/pictures/REVISION"
probis_th = 5
com_len_all = pd.DataFrame(columns=["Organism", "Cluster size", "# Cluster"])
singles_eval = []
results = []
pocket_th = 0.5
COM_DICT = {}
all_singles = []
pocket_df = []

for org in ["ECOLI", "YEAST", "CANAL", "ARATH", "ORYSJ", "MAIZE", "SOYBN", "DROME","CAEEL", "MOUSE", "HUMAN"]:
    bind_dict = pkl.load(open(f"../pm_clustering/results/bind_dict/{org}_pockets.pkl", "rb"))
    dom_dict = pkl.load(open(f"../pm_clustering/results/pocket_desc/{org}/dom_dict.pkl", "rb"))
    org_pockets = []
    for prot, pockets in bind_dict.items(): 
        for pocket in pockets:
            pocket_df.append(pocket)
            org_pockets.append(pocket)
    org_pockets = pd.DataFrame.from_records(org_pockets)
    probis_path = "../programs/probis"
    for pocket_th in [0.9, 0.8, 0.7, 0.6, 0.5]:
        print(f"Processing {org} with pocket threshold {pocket_th}")
        to_drop = org_pockets[org_pockets["prob"] < pocket_th]["pocket_id"].to_list()
        probis_df = pd.read_csv(f"{probis_path}/{org}/align_scores.csv", index_col = 0)

        probis_df.drop(to_drop, axis=0, inplace=True)
        probis_df.drop(to_drop, axis=1, inplace=True)

        align_scores = probis_df.copy()
        if probis_th > 1: 
            align_scores[align_scores == 1] = 0
        align_scores[align_scores >= probis_th] = 1
        align_scores[align_scores != 1] = 0

        n_components_probis, labels_probis = connected_components(csgraph = csr_matrix(align_scores), directed = False, return_labels = True)
        probis_count = Counter(labels_probis)
        n_probis_clust = len([k for k, v in probis_count.items() if v > 1])
        n_probis_singles = len([k for k, v in probis_count.items() if v == 1])

        num_bs = len(probis_df)

        singles = probis_df.loc[:, (probis_df == 0).all(axis=0)]
        singles = probis_df.loc[(probis_df == 0).all(axis=1)]
        singles = list(singles.index)
        num_singles = len(singles)

        probis_df = probis_df.loc[:, (probis_df != 0).any(axis=0)]
        probis_df = probis_df.loc[(probis_df != 0).any(axis=1)]
        # pocket_df["is_single"] = [True if pocket in singles.index else False for pocket in pocket_df["pocket_id"]]
        # probis_df = shuffle_symmetric_matrix(probis_df) 
        
        pocket_list = list(probis_df.index)
        edge_list = []
        for idx, pocket1 in enumerate(pocket_list): 
            for pocket2 in pocket_list[idx+1:]:
                score = probis_df.at[pocket1, pocket2]
                if score != 0 and not np.isnan(score):
                    edge_list.append(f"{pocket1} {pocket2} {score}\n")
        with open(f"edges/{org}_edges.txt", "w") as f:
            edge_list = "".join(edge_list)
            f.write(edge_list)

        g = nx.read_weighted_edgelist(open(f"edges/{org}_edges.txt", "r"))
        w_list = [edge[2]["weight"] for edge in g.edges(data=True)]
        coms = algorithms.leiden(g, weights=w_list)
        com_len = [len(com) for com in coms.communities]
        com_len = list(Counter(com_len).items())
        COM_DICT[org] = coms.communities

        # com_len.insert(0, (1, len(singles)))
        org_series = [org] * len(com_len)
        size, num = zip(*com_len)
        th_list = [pocket_th] * len(com_len)
        com_len = pd.DataFrame(list(zip(org_series, size, num, th_list)), columns = ["Organism", "Cluster size", "# Cluster", "Pocket_th"])
        com_len_all = pd.concat([com_len_all, com_len])
        # results.append([org, len(coms.communities), num_singles, round(num_singles/num_bs, 3), com_len["Cluster size"].max()])
        if pocket_th == 0.5:
            all_singles += singles

        print(f"Organism: {org}, #Communities: {len(coms.communities)}, #Singletons: {num_singles}, largest community: {com_len['Cluster size'].max()}")
        dom_df = pd.DataFrame(coms.communities[0], columns=["Pocket_ID"])
        dom_df["Domain"] = dom_df["Pocket_ID"].map(dom_dict)
        dom_df["Domain"] = [str(dom["domains"]) for dom in dom_df["Domain"]]
        dom_count = Counter(dom_df["Domain"])
        most_common = dom_count.most_common(2)
        num_kinase = 0 
        for dom, count in dom_count.items():
            if "Protein kinase domain" in dom: 
                num_kinase += count
        print(f"Most common protein domains in largest community: {most_common} with {round(num_kinase/len(dom_df), 3)} protein kinase domains")
        results.append([org, pocket_th, num_bs, len(coms.communities), num_singles, round(num_singles/num_bs, 3), com_len["Cluster size"].max(), most_common[0][0], most_common[0][1], most_common[1][0], most_common[1][1], round(num_kinase/len(dom_df), 3), n_probis_clust, n_probis_singles])
         

results = pd.DataFrame(results, columns=["Organism", "Pocket_th", "#Pockets", "#Communities", "#Singletons", "Single_frac", "Largest_comm", "Most_dom_1", "Num_1", "Most_dom_2", "Num_2", "Kinase_prec", "#Probis_cluster", "#Probis_singles"])
results.to_csv(f"{PLOT_DIR}/communities.csv", index=False)

# Evaluate singletons 
pocket_df = pd.DataFrame.from_records(pocket_df)
pocket_df["probis_single"] = ["Singletons" if pocket in all_singles else "Non-singletons" for pocket in pocket_df["pocket_id"]]
pocket_df["num_res"] = [len(res) for res in pocket_df["res"]]

singles_all = pocket_df[pocket_df["probis_single"] == "Singletons"]
probis_all = pocket_df[pocket_df["probis_single"] == "Non-singletons"] 

prop_dict = {"hydro": "Hydrophobicity", "aro": "Aromaticity", "net_charge": "Net charge",
            "sasa": "SASA [\u212B\u00b2]","plddt": "pLDDT", "pae": "PAE [\u212B]",
            "prob": "Probability", "num_res": "# Residues",
            }
titles = ["A", "B", "C", "D", "E", "F", "G", "H"]

dpi = 300
fig_width_px = 2250
fig_height_px = 2625
figsize_in = (fig_width_px / dpi, fig_height_px / dpi)
sns.set(style='white', font_scale=1)
n_rows = 4
n_cols = 2
row = 0
col = 0
fig, axes = plt.subplots(nrows=n_rows, ncols=n_cols, figsize=figsize_in, dpi=dpi)

ttest_results = []
for idx, (prop, label) in enumerate(prop_dict.items()):
    ax = axes[row, col]
    sns.violinplot(x='probis_single', y=prop, data=pocket_df, 
                    palette={'Non-singletons': 'lightseagreen', 'Singletons': 'darkorange'},
                    inner='box', linewidth=1.5, ax=ax)

    # Customize plot appearance
    ax.set_title(titles[idx], fontsize=12, weight='bold', loc='left')
    ax.set_ylabel(label, fontsize=10)
    ax.set_xlabel('')
    ax.tick_params("x", rotation=45)
    
    # Perform t-test 
    t_stat, p_val = ttest_ind(list(probis_all[prop]), list(singles_all[prop]))
    cohen_d, std_single, std_non_single = get_cohen_d(singles_all[prop], probis_all[prop])
    ttest_results.append([prop, t_stat, p_val, probis_all[prop].median(), singles_all[prop].median(), round(std_single,3), round(std_non_single,3), round(cohen_d,3)])

    col += 1
    if col == n_cols:
        col = 0
        row += 1

plt.tight_layout()
plt.savefig(f"{PLOT_DIR}/single_vs_non_single.tif",dpi=dpi, format="tif")

pocket_df.to_csv(f"{PLOT_DIR}/pocket_df.csv", index=False)
ttest_results = pd.DataFrame(ttest_results, columns=["Property", "t_stat", "p_val", "non_single_med", "single_med", "std_single", "std_non_single", "cohen_d"])
ttest_results.to_csv(f"{PLOT_DIR}/ttest_single_vs_non_singles.csv")

color_dict = {"ECOLI": "#0000FF", "CAEEL": "#40E0D0","YEAST": "gray", "CANAL": "#d62728", "DROME": "#9467bd",
             "MOUSE": "#f0e442", "HUMAN": "#FF9900", "ARATH": "#990099", "ORYSJ": "green",
             "SOYBN": "brown", "MAIZE": "black"}
 
com_len_all = com_len_all[com_len_all["Pocket_th"] == 0.5]
com_len_all["log(Cluster size)"] = np.log10(com_len_all["Cluster size"].to_list())
com_len_all["log(# Cluster)"] = np.log10(com_len_all["# Cluster"].to_list())

# Set plot style
sns.set(style="whitegrid", font_scale=1)

# Create the line plot
dpi = 300
fig_width = 4
fig_height = 4

fig = plt.figure(figsize=(fig_width, fig_height), dpi=dpi)
sns.lineplot(
    data=com_len_all,
    x="log(Cluster size)",
    y="log(# Cluster)",
    hue="Organism",
    palette=color_dict,
    marker="o",
)

# Customize appearance
plt.ylabel("log($N_{\mathrm{Clusters}}$)", fontsize=10)
plt.xlabel("log(Cluster size)", fontsize=10)
ax = plt.gca()
ax.set_facecolor("white")
plt.grid(False)
sns.despine()

# Adjust legend
ax.legend(title="Species")

# Save the figure
plt.tight_layout()
plt.savefig(f"{PLOT_DIR}/cluster_sizes.tif", dpi=300, format="tif")
pkl.dump(COM_DICT, open(f"{PLOT_DIR}/all_communities.pkl", "wb"))
