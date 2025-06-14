import os
import glob
import json
import pandas as pd
import numpy as np
import pickle as pkl 
from scipy.sparse.csgraph import connected_components
from scipy.sparse import csr_matrix


if __name__ == "__main__":
    
    POCKET_DIR = "../../pm_clustering/results/bind_dict"
    # ORGANISMS = ["ECOLI", "YEAST", "CANAL", "ARATH", "ORYSJ", "MAIZE", "DROME", "CAEEL", "MOUSE", "HUMAN"]
    ORGANISMS = ["SOYBN"]
    for org in ORGANISMS:
        print(org)
        srf_dir = f"{org}/results"
        out_file = f"{org}/align_scores.csv"
        pocket_dict = pkl.load(open(f"{POCKET_DIR}/{org}_pockets.pkl", "rb"))
        pocket_ids = [pocket["pocket_id"] for prot, pockets in pocket_dict.items() for pocket in pockets]
        align_scores = pd.DataFrame(np.zeros((len(pocket_ids), len(pocket_ids))), index=pocket_ids,columns=pocket_ids)
        clusters = pd.DataFrame(np.zeros((len(pocket_ids), len(pocket_ids))), index=pocket_ids,columns=pocket_ids)
        # z_scores = pd.DataFrame(index=list(id_dict.values),columns=list(id_dict.values))
        # e_values = pd.DataFrame(index=list(id_dict.values),columns=list(id_dict.values))
    
        for filename in glob.glob(f"{srf_dir}/*json"):
            pocket_id = os.path.basename(filename)[:-5]
            if pocket_id not in pocket_ids: 
                continue

            data = open(filename)
            data = json.load(data)

            # print(data[0]) 
            for align in data:
                pocket_id2 = align.get("pdb_id")
                if pocket_id2 not in pocket_ids: 
                    continue
                # print(srf_id2, pdb_id2)
                if align.get("pdb_id") == pocket_id:
                    clusters.at[pocket_id, pocket_id] = 0 

                else: 
                    align_score = align.get("alignment")[0].get("scores").get("alignment_score")
                    if align_score != None:
                        align_scores.at[pocket_id2, pocket_id] = align_score
                        align_scores.at[pocket_id, pocket_id2] = align_score
                        clusters.at[pocket_id2, pocket_id] = 1
                        clusters.at[pocket_id, pocket_id2] = 1

        tmp = csr_matrix(clusters.fillna(0))
        n_components, labels = connected_components(csgraph = tmp, directed = False, return_labels = True)
        print(f"Number of components: {n_components}")
        real_cluster = []
        max_count = None
        min_count = None
        for i in range(n_components):
            count = list(labels).count(i)
            if count > 1:
                real_cluster.append(count)
                if max_count == None or count > max_count[1]:
                    max_count = (i, count)
                if min_count == None or count < min_count[1]:
                    min_count = (i, count)
        # num_cluster.append(len(real_cluster))
        print(f"Number of real clusters: {len(real_cluster)}, Min length: {min_count}, Max length: {max_count}")

        clusters.to_csv(f"{org}/clusters.csv")
        align_scores.to_csv(out_file)
