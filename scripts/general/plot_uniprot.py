import os
import glob
import yaml
import pandas as pd 
import pickle as pkl

import matplotlib.pyplot as plt
import seaborn as sns

from scipy.stats import fisher_exact
from statsmodels.stats.multitest import fdrcorrection


if __name__ == "__main__":
    # Load config and define parameters
    with open("../../config.yaml", "r") as f:
        config = yaml.safe_load(f)

    PLOT_DIR = config["paths"]["plots"]
    RESULTS_PATH = config["paths"]["results"]
    STRUC_PATH = config["paths"]["proteins"]
    BIND_DICT_PATH = config["paths"]["bind_dict"]
    KNOWN_POCKETS_PATH = config["paths"]["known_pockets"]
    DB_PATH = config["paths"]["db_files"]
    SEQ_DICT = pkl.load(open(f"{STRUC_PATH}/seq_len_dict.pkl", "rb"))
    FS_PATH = f"{STRUC_PATH}/fs_clust_proteins.pkl"
    FS_DICT = pkl.load(open(FS_PATH, "rb"))

    SPECIES = [species["id"] for species in config["species"]]
    class_color = config["ligand_color"]
    
    class_map = {"CHEBI:16646": "Carbohydrate & derivatives", "CHEBI:63299": "Carbohydrate & derivatives", "CHEBI:167559": "Glycan", "CHEBI:37163": "Glycan",
                 "CHEBI:35381": "Monosaccharide & derivatives", "CHEBI:63367": "Monosaccharide & derivatives", "CHEBI:15841": "Oligo- & polypeptide",
                 "CHEBI:25676": "Oligo- & polypeptide", "CHEBI:33709": "AA & derivatives", "CHEBI:35238": "AA & derivatives",
                "CHEBI:83821": "AA & derivatives", "CHEBI:16991": "Nucleic acid", "CHEBI:33697": "Nucleic acid", "CHEBI:18282": "Nucleobase/nucleoside/nucleotide",
                "CHEBI:33838": "Nucleobase/nucleoside/nucleotide", "CHEBI:36976": "Nucleobase/nucleoside/nucleotide", "CHEBI:18059": "Lipid",
                "CHEBI:36914": "Inorganic ion", "CHEBI:24867": "Inorganic ion", "CHEBI:30413": "Heme", "CHEBI:33733": "Hetero nuclear cluster",
                "Other compound": "Other compound"}
    
    chebi_dict = {lig: [] for lig in class_color.keys()}
    for chebi, lig in class_map.items():
        chebi_dict[lig].append(chebi)
    
    # Load known pockets and predicted pockets
    PRED_POCKETS = []
    for species in SPECIES:
        pred_dict = pkl.load(open(f"{BIND_DICT_PATH}/{species}_pockets.pkl", "rb"))
        for prot, pockets in pred_dict.items():
            for pocket in pockets:
                pocket["prot"] = prot
                PRED_POCKETS.append(pocket)
    PRED_POCKETS = pd.DataFrame.from_records(PRED_POCKETS)
    PRED_POCKETS["#res"] = PRED_POCKETS["res"].apply(lambda x: len(x))
    PRED_POCKETS["#res"] = PRED_POCKETS["#res"].astype(int)

    known_dict = pkl.load(open(f"{KNOWN_POCKETS_PATH}/known_pockets_up.pkl", "rb"))
    KNOWN_POCKETS = []
    for prot, pockets in known_dict.items():
        for pocket in pockets: 
            pocket["prot"] = prot 
            KNOWN_POCKETS.append(pocket)
    KNOWN_POCKETS = pd.DataFrame.from_records(KNOWN_POCKETS)
    KNOWN_POCKETS["found"] = [True if prob >= 0.5 else False for prob in KNOWN_POCKETS["prob"] ]
    KNOWN_POCKETS["#res"] = KNOWN_POCKETS["res"].apply(lambda x: len(x))
    KNOWN_POCKETS["#res"] = KNOWN_POCKETS["#res"].astype(int)

    pocket_info = []
    for species in SPECIES:
        up_path = os.path.join(DB_PATH, species)
        struc_path = os.path.join(STRUC_PATH, species)

        tm_dict = pkl.load(open(f"{RESULTS_PATH}/pocket_desc/{species}/{species}_tm_info.pkl", "rb"))
        tm_failed = [prot for prot, tests in tm_dict.items() if tests["num_passed"] < 6]
        if os.path.isfile(os.path.join(up_path, f"frag_df.tsv.gz")):
            fragments = list(pd.read_csv(os.path.join(up_path, f"frag_df.tsv.gz"), index_col="Entry", delimiter="\t", compression="infer").index)
        else: 
            fragments = []
        all_prots = [os.path.basename(prot).split("-")[1] for prot in glob.glob(f"{struc_path}/*pdb")]
        all_prots = [prot for prot in all_prots if prot not in fragments and SEQ_DICT.get(prot, 0) >= 100]
        all_prots = [prot for prot in all_prots if prot not in tm_failed]
        fs_cluster = [FS_DICT.get(prot)["clusterID"] for prot in all_prots if FS_DICT.get(prot) is not None]
        num_fs = len(set(fs_cluster)) if fs_cluster else 0
        num_prot = len(all_prots)

        tmp_pred = PRED_POCKETS[PRED_POCKETS["species"] == species]
        num_pred = len(tmp_pred)
        tmp_known = KNOWN_POCKETS[KNOWN_POCKETS["species"] == species]
        num_known = len(tmp_known)

        # ECO in the dataset
        # ECO:0000250 - sequence similarity evidence used in manual assertion
        # ECO:0000255 - match to sequence model evidence used in manual assertion
        # ECO:0000256 - match to sequence model evidence used in automatic assertion

        # ECO:0000269 - experimental evidence used in manual assertion
        # ECO:0000305 - curator inference used in manual assertion
        # ECO:0007744 - combinatorial computational and experimental evidence used in manual assertion
        # ECO:0007829 - combinatorial computational and experimental evidence used in automatic assertion

        # Get all known pockets with experimental evidence
        num_exp = len(KNOWN_POCKETS[(KNOWN_POCKETS["species"] == species) & (KNOWN_POCKETS["evidence"] == "ECO:0000269")])
        # Get all known pockets with no sequence similarity evidence
        no_seq_sim = len(KNOWN_POCKETS[(KNOWN_POCKETS["species"] == species) & ~(KNOWN_POCKETS["evidence"].isin(["ECO:0000250", "ECO:0000255", "ECO:0000256"]))])
        seq_sim = len(KNOWN_POCKETS[(KNOWN_POCKETS["species"] == species) & (KNOWN_POCKETS["evidence"].isin(["ECO:0000250", "ECO:0000255", "ECO:0000256"]))])
        pocket_info.append({"species": species, "num_prot": num_prot, "num_pred": num_pred,
                            "num_fs": num_fs, "num_known": num_known,"num_exp": num_exp, 
                            "no_seq_sim": no_seq_sim, "seq_sim": seq_sim})
    pocket_info = pd.DataFrame.from_records(pocket_info)
    pocket_info.to_csv(f"{RESULTS_PATH}/known_pockets_info.csv", index=False)

    # Evaluate enrichment of found pockets per ligand class
    lig_classes = list(class_color.keys())
    pval_dict = {}
    results_dict = {}
    for lig in lig_classes: 
        num_pockets = len(KNOWN_POCKETS[KNOWN_POCKETS["class_name"] == lig])
        if num_pockets == 0:
            continue
        lig_found = len(KNOWN_POCKETS[(KNOWN_POCKETS["class_name"] == lig) & (KNOWN_POCKETS["found"] == True)])
        other_found = len(KNOWN_POCKETS[(KNOWN_POCKETS["class_name"] != lig) & (KNOWN_POCKETS["found"] == True)])

        lig_not_found = len(KNOWN_POCKETS[(KNOWN_POCKETS["class_name"] == lig) & (KNOWN_POCKETS["found"] == False)])
        other_not_found = len(KNOWN_POCKETS[(KNOWN_POCKETS["class_name"] != lig) & (KNOWN_POCKETS["found"] == False)])

        conf_matrix = [[lig_found, other_found], [lig_not_found, other_not_found]]
        fold, pvalue = fisher_exact(conf_matrix)
        
        pval_dict[lig] = pvalue
        results_dict[lig] = {"chebi_ids": chebi_dict.get(lig, []), "n_pockets": num_pockets, "found": lig_found, "not_found": lig_not_found, "fold": round(fold, 2), "pval": pvalue}

    rejected, padj = fdrcorrection(list(pval_dict.values()), method="n")
    for idx, lig in enumerate(pval_dict.keys()):
        results_dict[lig]["padj"] = padj[idx]
        results_dict[lig]["rejected"] = rejected[idx]
    found_enrichment = pd.DataFrame.from_dict(results_dict, orient="index")
    found_enrichment = found_enrichment.sort_values(by="found", ascending=False)
    found_enrichment.to_csv(f"{PLOT_DIR}/ligands_found_enrichment.csv")

    # Plot pocket size distribution
    size_df = pd.DataFrame({
        "Pocket type": ["Known"] * len(KNOWN_POCKETS),
        "# residues": KNOWN_POCKETS["#res"]
    })
    tmp = pd.DataFrame({
        "Pocket type": ["Predicted"] * len(PRED_POCKETS),
        "# residues": PRED_POCKETS["#res"]
    })
    size_df = pd.concat([size_df, tmp])

    dpi = 300
    fig_width = 3.3
    fig_height = 3.3

    plt.figure(figsize=(fig_width, fig_height), dpi=dpi)
    sns.violinplot(
        x="Pocket type",
        y="# residues",
        data=size_df,
        inner="box",
        hue="Pocket type",
        palette={"Known": "lightseagreen", "Predicted": "darkorange"},
        linewidth=1
    )

    plt.xlabel("Pocket type", fontsize=10)
    plt.ylabel("# residues", fontsize=10)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.tight_layout()
    plt.savefig(f"{PLOT_DIR}/pocket_size.tif", dpi=dpi, format='tif', bbox_inches='tight')
    plt.close()
