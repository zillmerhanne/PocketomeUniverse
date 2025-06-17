import os
import glob
import json
import yaml
import pandas as pd
import numpy as np
import pickle as pkl 
from collections import Counter
from scipy.sparse.csgraph import connected_components
from scipy.sparse import csr_matrix


if __name__ == "__main__":
    # Load configurations
    with open("config.yaml", "r") as f:
        config = yaml.safe_load(f)
    BIND_DICT_PATH = config["paths"]["bind_dict"]
    SPECIES = [species["id"] for species in config["species"]]
    PROBIS_PATH = config["paths"]["probis_results"]
    PROBIS_TH = config["general"]["probis_th"]

    RESULTS = []
    for species in SPECIES:
        print("Processing species:", species)
        srf_dir = f"{PROBIS_PATH}/{species}/srf"
        out_file = f"{PROBIS_PATH}/{species}/align_scores.csv.gz"
        bind_dict = pkl.load(open(f"{BIND_DICT_PATH}/{species}_pockets.pkl", "rb"))
        pocket_ids = [pocket["pocket_id"] for prot, pockets in bind_dict.items() for pocket in pockets]
        align_scores = pd.DataFrame(np.zeros((len(pocket_ids), len(pocket_ids))), index=pocket_ids,columns=pocket_ids)
        # z_scores = pd.DataFrame(index=list(id_dict.values),columns=list(id_dict.values))
        # e_values = pd.DataFrame(index=list(id_dict.values),columns=list(id_dict.values))

        # Read all json files in the directory
        # Extract aligned pockets and their scores
        for filename in glob.glob(f"{srf_dir}/*json"):
            pocket_id = os.path.basename(filename)[:-5]
            if pocket_id not in pocket_ids: 
                continue
            with open(filename) as data_file:
                data = json.load(data_file)
            for align in data:
                pocket_id2 = align.get("pdb_id")
                if pocket_id2 not in pocket_ids: 
                    continue
                
                if align.get("pdb_id") != pocket_id: 
                    align_score = align.get("alignment")[0].get("scores").get("alignment_score")
                    if align_score != None:
                        align_scores.at[pocket_id2, pocket_id] = align_score
                        align_scores.at[pocket_id, pocket_id2] = align_score

        # Simple clustering based on alignment scores 
        # Convert alignment scores to binary clusters based on threshold
        # Determine connected components in the binary matrix
        clusters = align_scores.copy()
        clusters[clusters < PROBIS_TH] = 0
        clusters[clusters >= PROBIS_TH] = 1
        tmp = csr_matrix(clusters)
        n_components, labels = connected_components(csgraph = tmp, directed = False, return_labels = True)
        print(f"Number of components: {n_components}")
        cluster_count = Counter(labels)
        real_cluster = len([count for count in cluster_count.values() if count > 1])
        singletons = [count for count in cluster_count.values() if count == 1]
        max_count = cluster_count.most_common(1)[0][1]
        print(f"Number of real clusters: {len(real_cluster)}, Max length: {max_count}")
        RESULTS.append({
            "species": species,
            "n_components": n_components,
            "real_clusters": real_cluster,
            "largest_cluster": max_count,
            "singletons": len(singletons)
        })
        # Save results 
        align_scores.to_csv(out_file)

    # Save summary results
    summary_df = pd.DataFrame.from_records(RESULTS)
    summary_file = f"{PROBIS_PATH}/probis_clustering_th={PROBIS_TH}.csv.gz"
    summary_df.to_csv(summary_file, index=False)