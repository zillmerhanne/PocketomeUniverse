import os
import numpy as np
import pandas as pd
import pickle as pkl
import plotly.express as px
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA, FastICA
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import DBSCAN
from sklearn.metrics import silhouette_score, pairwise_distances
from sklearn.manifold import TSNE
from scipy.stats import pearsonr, spearmanr
import umap 
from scipy.spatial.distance import pdist
import matplotlib.pyplot as plt 
import seaborn as sns


def plot_dim_red(df, red_type, scatter_configs, out_file):
    sns.set(style="white")
    dpi = 300
    fig_width_px = 2050
    fig_height_px = 2625
    figsize_in = (fig_width_px / dpi, fig_height_px / dpi)
    fig, axes = plt.subplots(2, 2, figsize=figsize_in, dpi=dpi)
    plt.subplots_adjust(hspace=0.3, wspace=0.3)
    
    # Shared axis labels
    x_col = f"{red_type}1"
    y_col = f"{red_type}2"
    x_pos = [0, 0, 1, 1]
    y_pos = [0, 1, 0, 1]
    # Define the plots
    for i, (color_col, plot_label, cbar_title, cmap) in enumerate(scatter_configs):
        sc = sns.scatterplot(
            x=x_col, y=y_col,
            hue=color_col,
            palette=cmap,
            data=df,
            ax=axes[x_pos[i], y_pos[i]],
            s=2,
            legend=False
        )
        axes[x_pos[i], y_pos[i]].set_title(plot_label, fontsize=16, weight='bold', loc='left')
        axes[x_pos[i], y_pos[i]].set_xlabel(x_col, fontsize=10)
        axes[x_pos[i], y_pos[i]].set_ylabel(y_col, fontsize=10)
    
        # Add colorbar manually
        norm = plt.Normalize(df[color_col].min(), df[color_col].max())
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=axes[x_pos[i], y_pos[i]], orientation='vertical', fraction=0.046, pad=0.04)
        cbar.ax.tick_params(labelsize=10)
        cbar.set_label(cbar_title, fontsize=10)
    
    # Save the figure
    fig.tight_layout()
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.savefig(out_file, dpi=dpi, bbox_inches='tight')
    plt.close()
    return 0

def get_embeddings(results_path):
    results = pkl.load(open(results_path, "rb"))
    pocket_embeddings = {}
    for pocket1, pocket2 in results["pairs"]:
        pocket_id1 = os.path.basename(pocket1["pocket"])[:-4]
        pocket_id2 = os.path.basename(pocket2["pocket"])[:-4]
        if pocket_id1 not in pocket_embeddings and pocket_id1 in pocket_ids:
            pocket_embeddings[pocket_id1] = pocket1["descriptor"].numpy()[0]
        if pocket_id2 not in pocket_embeddings and pocket_id2 in pocket_ids:
            pocket_embeddings[pocket_id2] = pocket2["descriptor"].numpy()[0]
        
    return pocket_embeddings


def k_distance_plot(data,th):
  """Calculates and plots the k-distance graph.

  Args:
    data: The data points.
    k: The number of nearest neighbors.
  """
  k_list = [3, 5, 7, 9, 11]
  fig, ax = plt.subplots()
  for k in k_list:
    neigh = NearestNeighbors(n_neighbors=k+1)
    nbrs = neigh.fit(data)
    distances, indices = nbrs.kneighbors(data)
    distances = np.sort(distances[:, k], axis=0)
  
    ax.plot(distances, label=k)
  plt.xlabel('Data Points')
  plt.ylabel('Distance')
  plt.title(f'K-Distance Graph')
  plt.legend()
  plt.savefig(f"k_dist_{th}.png")


def cluster_DBSCAN(org, th, embeddings_np):
    results = []
    for eps in [0.05, 0.1, 0.15, 0.2]:
        for min_samples in [3, 5, 7, 9]:
            dbscan = DBSCAN(eps=eps, min_samples=min_samples)  # Adjust eps and min_samples as needed
            labels = dbscan.fit_predict(embeddings_np)

            # Evaluate cluster quality
            silhouette_avg = silhouette_score(embeddings_np, labels)
            print("Silhouette coefficient:", silhouette_avg)
            # Analyze cluster results
            n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
            unique, counts = np.unique(labels, return_counts = True)
            counts = dict(zip(unique, counts))
            print(f"#Clusters:in {n_clusters}, #Singletons: {counts.get(-1)}")

            results.append([org, th, eps, min_samples, silhouette_avg, n_clusters, counts.get(-1)])
    return results 


def get_correl(embed, desc_np, desc_name, corr_type):
    results = [desc_name, corr_type]
    for idx in range(embed.shape[1]):
        vect = embed[:, idx]
        # print(len(vect))
        if corr_type == "Pearson":
            corr, p = pearsonr(desc_np, vect)
        else:
            corr, p = spearmanr(desc_np, vect)   
        results.append(corr)
      
    for feat in ["PC1", "PC2", "PC3", "PC4","IC1", "IC2", "IC3"]:
        vect = pockets_all[feat].to_numpy()
        if corr_type == "Pearson": 
            corr, p = pearsonr(desc_np, vect)
        else: 
           corr, p = spearmanr(desc_np, vect)
        results.append(corr)
    return results


if __name__ == "__main__":
    cluster_df = []
    # Add the definitions to a config file 
    results_path = ".."
    organisms = ["ECOLI", "CANAL", "YEAST", "HUMAN", "MOUSE", "ARATH", "DROME", "CAEEL", "MAIZE", "ORYSJ", "SOYBN"]
    DEEPLYTOUGH_PATH = ".."
    RESULTS_PATH = "../../pm_clustering/results"
    PLOT_DIR = f"{RESULTS_PATH}/pictures/REVISION"
    DB_PATH = "../../../Pocketeome/db_files"
    SEQ_PATH = "../../pm_clustering/results/proteins/seq_len_dict.pkl"
    SEQ_DICT = pkl.load(open(SEQ_PATH, "rb"))
    scatter_configs = [
        ("hydro", "A","Hydrophobicity", "RdBu"),
        ("aro", "B", "Aromaticity", "viridis"),
        ("net_charge", "C","Net charge", "RdBu")
        ("sasa", "D", "SASA [\u212B\u00b2]", "viridis"),
        ("plddt", "E", "pLDDT", "viridis"),
        ("pae", "F", "PAE [\u212B]", "viridis")
        ("#res", "G", "# residues", "viridis")
        ("prob", "H", "Probability", "viridis"),
    ]

    
    with open("../../pm_clustering/data/fasta/clustered_proteins.txt", "r") as f: 
        lines = f.readlines()
        clustered = [line.strip() for line in lines[1:]]

    th = 0.5 

    embeddings_all = {}
    pockets_all = []
    pocket_ids = []
    for org in organisms:   
        print(f"Load data for: {org}")
        bind_dict = pkl.load(open(f"../../pm_clustering/results/bind_dict/{org}_pockets.pkl", "rb"))
        tm_dict = pkl.load(open(f"../../pm_clustering/results/pocket_desc/{org}/{org}_tm_info.pkl", "rb"))
        
        pockets_all += [pocket for pockets in bind_dict.values() for pocket in pockets]
        pocket_ids += [pocket["pocket_id"] for pockets in bind_dict.values() for pocket in pockets]
        deeplytough_result = f"{DEEPLYTOUGH_PATH}/{org}/Custom-DeeplyTough-networks.pickle"
        embeddings = get_embeddings(deeplytough_result)
        embeddings_all.update(embeddings)
        
    
    pockets_all = pd.DataFrame.from_records(pockets_all)
    pockets_all["#res"] = pockets_all["res"].apply(lambda x: len(x))
    pockets_all.sort_values(by="pocket_id", inplace=True)
    embeddings_all = dict(sorted(embeddings_all.items(), key=lambda item: item[0]))
    pkl.dump(embeddings_all, open("all_embeddings.pkl", "wb"))
    # org_map_all = dict(sorted(org_map_all.items(), key=lambda item: item[0])) 
    
    embeddings_np = np.array(list(embeddings_all.values())) 
    # dist = pairwise_distances(embeddings_np, n_jobs=10)
    # with open("pairwise_dist.npy", "wb") as f: 
    #    np.save(f, dist)              
    
    scaler = StandardScaler()
    embeddings_norm = scaler.fit_transform(embeddings_np)

    print("Calculate PCA")
    pca = PCA(n_components=4)
    pca_out = pca.fit_transform(embeddings_np)
    print("Explained variance for normalized PCA results")
    for idx, explained_var in enumerate(pca.explained_variance_ratio_[:5]):
        print(f"PC{idx+1}: {round(explained_var,2)}")
            
    pockets_all["PC1"] = pca_out[:, 0]
    pockets_all["PC2"] = pca_out[:, 1]
    pockets_all["PC3"] = pca_out[:, 2]
    pockets_all["PC4"] = pca_out[:, 3]
    
    plot_dim_red(pockets_all, "PC", scatter_configs, f"{PLOT_DIR}/pca_desc.tif")

    n_neighbors = [10, 20, 50, 100, 200] 
    for n in n_neighbors:       
        print("Calculate UMAP")
        reducer = umap.UMAP(
                n_neighbors = n, # local vs global clustering, with low n_neighbors the focus is more on local shapes (~2)
                min_dist=0,       # low because we are interested in clustering 
                n_components=2,   # to plot in 2D
                random_state=1
                )
        umap_out = reducer.fit_transform(embeddings_np)
        pockets_all[f"UMAP1_n{n}"] = umap_out[:,0]
        pockets_all[f"UMAP2_n{n}"] = umap_out[:,1]

        pockets_all[f"UMAP1"] = umap_out[:,0]
        pockets_all[f"UMAP2"] = umap_out[:,1]
        plot_dim_red(pockets_all, "UMAP", scatter_configs, f"{PLOT_DIR}/umap_n{n}_desc.tif")
        
    pockets_all.drop(columns=["UMAP1", "UMAP2"], inplace=True)
            
    print("Calculate tSNE")
    perplexities = [10, 20, 30, 40, 50, 100]
    for perplexity in perplexities:
        # see: https://scikit-learn.org/stable/auto_examples/manifold/plot_t_sne_perplexity.html for effect of perplexity
        tsne = TSNE(n_components=2, 
                    learning_rate="auto", 
                    init='random', 
                    perplexity=perplexity, 
                    n_iter=5000, 
                    random_state=1)
        tsne_out = tsne.fit_transform(embeddings_np)
        print(f"Number of iterations needed for tSNE run: {tsne.n_iter_}")
        pockets_all[f"tSNE1_p{perplexity}"] = tsne_out[:,0]
        pockets_all[f"tSNE2_p{perplexity}"] = tsne_out[:,1]
#
        pockets_all[f"tSNE1"] = tsne_out[:,0]
        pockets_all[f"tSNE2"] = tsne_out[:,1]
        plot_dim_red(pockets_all, "tSNE", scatter_configs, f"{PLOT_DIR}/tsne_p{perplexity}_desc.tif")
        
    pockets_all.drop(columns=["tSNE1","tSNE2"], inplace=True)

    # print("Calculate tSNE")
    # tsne = TSNE(n_components=2, learning_rate="auto", init='random', perplexity=15, n_iter=5000, random_state=1)   # see: https://scikit-learn.org/stable/auto_examples/manifold/plot_t_sne_perplexity.html for effect of perplexity
    
    # tsne_out = tsne.fit_transform(np.array([pca_out.T[0], pca_out.T[1], pca_out.T[2], pca_out.T[3], pca_out.T[4]]).T)
    # tsne_out = pd.DataFrame(list(zip(list(org_map_all.values()),list(embeddings_all.keys()),tsne_out.T[0], tsne_out.T[1])), columns=["Organism", "Pocket ID","tSNE1", "tSNE2"])
    # tsne_out.to_csv("tsne_on_pca.csv", index=False)
    # print(f"Number of iterations needed for tSNE run: {tsne.n_iter_}")
    # dim_red_df["tSNE1"] = tsne_out[:,0]
    # dim_red_df["tSNE2"] = tsne_out[:,1]
    
    ica = FastICA(n_components=3)
    ica_out = ica.fit_transform(embeddings_np)
 
    pockets_all["IC1"] = ica_out[:,0]
    pockets_all["IC2"] = ica_out[:,1]
    pockets_all["IC3"] = ica_out[:,2]
    plot_dim_red(pockets_all, "IC", scatter_configs, f"{PLOT_DIR}/ica_desc.tif")
    pockets_all.to_csv(f"dim_red_all.csv", index=False)

    columns = ["Pocket Desc", "Type"]
    columns += [f"Feat{idx+1}" for idx in range(embeddings_np.shape[1])]
    columns += ["PC1", "PC2", "PC3", "PC4", "IC1", "IC2", "IC3"]
    
    feats = ["aro", "hydro", "plddt", "sasa", "prob", "net_charge", "pae", "#res"]
    correl_types = ["Pearson", "Spearman"]
    
    results_all = []
    for feat in feats: 
        for correl in correl_types: 
            results_all.append(get_correl(embeddings_np, pockets_all[feat].to_numpy(), feat, correl))
    results_all = pd.DataFrame(results_all, columns = columns)
    
    for dim_red in ["PC1", "PC2", "PC3", "PC4", "IC1", "IC2", "IC3"]:
        print(dim_red)
        for pocket_desc, corr_type, corr in zip(results_all["Pocket Desc"], results_all["Type"], results_all[dim_red]):
            if corr_type == "Spearman":
                print(pocket_desc, round(corr, 3))
    results_all.to_csv("Correlations.csv", index=False)
