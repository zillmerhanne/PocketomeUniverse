#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd 
import numpy as np 
import os 
import glob
import pickle
import yaml
from concurrent.futures import ProcessPoolExecutor, as_completed

from Bio.SeqUtils.ProtParam import ProteinAnalysis 
from Bio.SeqUtils import seq1
from Bio.PDB import PDBParser, Select, PDBIO

import seaborn as sns 
import matplotlib.pyplot as plt

from scipy.stats import kruskal


def get_pred_binding_sites(in_path, prot, th_pocket=0.3, greater_than=4, smaller_than=1000):
    """
    Get predicted binding sites from prediction files.

    Parameters:
        in_path (str): Path to the input files.
        prot (str): Protein identifier.
        th_pocket (float): Threshold for pocket probability.
        greater_than (int): Minimum number of residues in pocket.
        smaller_than (int): Maximum number of residues in pocket.

    Returns:
        pockets (list): List of pockets with residues and chains.
        probabilities (list): List of pocket probabilities.
    """
    path_options = [
        f"{in_path}/{prot}.pdb_predictions.csv",
        f"{in_path}/{prot}-F1-model_v1.pdb.gz_predictions.csv",
        f"{in_path}/{prot}-F1-model_v4.pdb.gz_predictions.csv"
    ]
    path = next((p for p in path_options if os.path.isfile(p)), None)
    if not path:
        return [], []

    pockets_df = pd.read_csv(path)
    pockets = []
    chains = []
    probabilities = []

    for i, pocket in enumerate(pockets_df[" residue_ids"]):
        prob = pockets_df.iloc[i, 3]
        if prob < th_pocket:
            break
        pocket_res = []
        chain = []
        for res in pocket.split():
            res_chain, res_id = res.split("_")
            res_id = int(res_id)
            chain.append(res_chain)
            pocket_res.append(res_id)
        
        if greater_than < len(pocket_res) < smaller_than:
            pockets.append(pocket_res)
            chains.append(chain)
            probabilities.append(prob)
        # print(pocket_res)
    if not pockets:
        return [], []

    return [list(zip(pocket, chain)) for pocket, chain in zip(pockets, chains)], probabilities


class ResidueSelector(Select):
    def __init__(self, residues_to_select):
        self.residues_to_select = residues_to_select
  
    def accept_residue(self, residue):
        if (residue.get_id()[1], residue.parent.id) in self.residues_to_select:
            return 1
        else:
            return 0


def create_pocket_pdb(pocket_path, prot, struc_file, pocket, num):
    """
    Create a PDB file for the pocket residues.
    Parameters:
        pocket_path (str): Path to save the pocket PDB file.
        prot (str): Protein identifier.
        struc_file (str): Path to the original PDB file.
        pocket (list): List of tuples containing residue IDs and chain IDs.
        num (int): Pocket number.
    Returns:
        int: Status code (0 for success).
    """
    out_path = os.path.join(pocket_path, f"{prot}_{num}.pdb")
    # Initialise a parser
    parser = PDBParser(QUIET=True)
    # Read the structure
    structure = parser.get_structure(prot, struc_file)
    # This assumes you have only one model in the PDB file
    io = PDBIO()
    io.set_structure(structure)
    io.save(out_path, ResidueSelector(pocket))
    return 0


def process_file(species, file, pred_path, th, out_path):
    """
    Process a single PDB file to extract binding sites and pLDDT values for each protein.
    Parameters:
        species (str): Species name.
        file (str): Path to the PDB file.
        pred_path (str): Path to the directory containing prediction files.
        th (float): Threshold for pocket probability.
        out_path (str): Path to save the output files.
    Returns:
        str: Protein identifier.
        dict: Dictionary containing the pockets and their properties.
        float: Mean pLDDT value for the protein.
    """
    result = {
        'has_pockets': 0,
        'pockets': [],
        'num_pockets': 0,
        "is_invalid": False
    }
    file = os.path.basename(file) 
    prot = "-".join(file.split("-")[:2]) 
    prot_pockets, prob = get_pred_binding_sites(pred_path, prot, th_pocket=th, greater_than=4, smaller_than=1000)
    plddt_dict, seq_dict = get_plddt(struc_path, prot)
    tm_info = tm_dict.get(prot[3:], None)
    plddt_mean = np.mean(list(plddt_dict.values()))
    pae_file = os.path.join(pae_path, f"{prot[3:]}_pae.npz")

    if len(plddt_dict) < 100 or not os.path.isfile(pae_file) or prot[3:] in fragments:
            result['is_invalid'] = True
            return prot, result, None
    if tm_info and tm_info["num_passed"] < 6:
            result['is_invalid'] = True
            return prot, result, None 

    if prot_pockets:        
        pae = get_pae(pae_file)
        result['has_pockets'] = 1
        for pocket_num, pocket in enumerate(prot_pockets):
            pocket_res, _ = zip(*pocket)
            try:
                pocket_plddt = np.mean([plddt_dict[res] for res in pocket_res])
                pocket_pae  = np.mean([pae[res1-1, res2-1] for res1 in pocket_res for res2 in pocket_res if res1 != res2])
            
                pocket_seq = [seq_dict[res] for res in pocket_res]
            except: 
                print(prot, pocket_res)
                return prot, result, None
            hydro, aro, net = get_hydro_aro(pocket_seq)
            create_pocket_pdb(out_path, prot, file, pocket, pocket_num)
            result['pockets'].append({"pocket_id": f"{prot}_{pocket_num}",
                                      "species": species,
                                      "res": pocket_res, "seq": pocket_seq,
                                      "prob": prob[pocket_num],
                                      "plddt": pocket_plddt,
                                      "pae": pocket_pae,
                                      "hydro": hydro,
                                      "aro": aro,
                                      "net_charge": net})
            result['num_pockets'] += 1                                                                                                 
    return prot, result, plddt_mean

def process_files_in_parallel(struc_path, pred_path, th, out_path, species, max_workers=None):
    """
    Process multiple PDB files in parallel to extract binding sites and pLDDT values.
    Parameters:
        struc_path (str): Path to the directory containing PDB files.
        pred_path (str): Path to the directory containing prediction files.
        th (float): Threshold for pocket probability.
        out_path (str): Path to save the output files.
        species (str): Species name.
        max_workers (int): Maximum number of workers for parallel processing.
        Returns:
            dict: Dictionary containing the pockets for each protein.
            int: Number of processed and valid proteins.
            dict: Dictionary containing pLDDT values for each protein.
    """
    files = glob.glob(os.path.join(struc_path, "*.pdb"))
    results = {}
    plddt_all = {}

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(process_file, species, file, pred_path, th, out_path): file for file in files}
        num_prot = 0                                                                                                                      
        for future in as_completed(futures):
            file = futures[future]
            # try:
            prot, result, plddt_mean = future.result()
            if result['is_invalid']:
                continue

            plddt_all[prot] = {"species": species,
                                "prot_plddt": plddt_mean, 
                                "has_pockets": result['has_pockets'],
                                }
            num_prot += 1
            if result['has_pockets']:
                results[prot] = result['pockets']
            # except Exception as exc:
            #     print(f'{file} generated an exception: {exc}')
                                                                                                                                                                                                                                                                                                    
        return results, num_prot, plddt_all
                                                                                                                                                                                                                                                                                                
def get_plddt(in_path, prot_id): 
    """
    Get the pLDDT values which are stored in the b factor position.

    Parameters:
        in_path (str): Path to the input files.
        prot_id (str): Protein identifier.

    Returns:
        str: Main chain identifier.
    """
    path = os.path.join(in_path, f"{prot_id}.pdb")
    plddt = {}
    seq = {}
    with open(path, "r") as f:
        for line in f:
            if line.startswith("ATOM") and line[12:15].strip() == "CA":
                plddt[int(line[22:26].strip())] = float(line[60:66].strip())
                seq[int(line[22:26].strip())] = line[17:20].strip()

    return plddt, seq


def get_pae(pae_path):  
    """
    Load the PAE data from a .npz file.
    Parameters:
        pae_path (str): Path to the PAE file.
    Returns:
        np.ndarray: PAE data array.
    """
    return np.load(pae_path)["data_array"]
    
def get_net_charge(seq):
    """
    Calculate the net charge of a protein sequence based on the amino acid composition.
    Parameters:
        seq (str): Protein sequence.
    Returns:
        net_charge (float): Net charge of the sequence.
    """
    # Net charge values for each amino acid
    # Special case: Histidine is set to 0.5
    NET_DICT = {"A": 0, "R": 1.0, "N": 0, "D": -1.0, "C": 0, "E": -1.0,
                "Q": 0, "G": 0, "H": 0.5, "I": 0, "L": 0, "K": 1.0,
                "M": 0, "F": 0, "P": 0, "S": 0, "T": 0, "W": 0,
                "Y": 0, "V": 0}
    net_charge = 0
    for aa in seq:
        net_charge += NET_DICT.get(aa)
    return net_charge

# Get the hydrophobicity and aromaticity of a pocket sequence
def get_hydro_aro(seq):
    """
    Calculate the hydrophobicity and aromaticity of a protein sequence.
    Parameters:
        seq (str): Protein sequence.
    Returns:
        hydro (float): Hydrophobicity value.
        aro (float): Aromaticity value.
        net (float): Net charge value.
    """
    seq = seq1("".join(seq))
    analysis = ProteinAnalysis(seq)
    hydro = round(analysis.gravy(), 3)
    aro = round(analysis.aromaticity(), 3)
    net = get_net_charge(seq)
    return hydro, aro, net


if __name__ == "__main__":

    # Load configuration file and define parameters
    with open("../../config/config.yaml", "r") as f:
        config = yaml.safe_load(f)

    TH = config["generalt"]["pocket_th"]
    THREADS = config["generalt"]["n_threads"]
    PAE_TH = config["generalt"]["pae_th"]
    PLDDT_TH = config["generalt"]["plddt_th"]

    SPECIES = [species["id"] for species in config["species"]]
    SPECIES_TITLES = [f"<b>{species}</b>" for species in SPECIES]

    STRUC_PATH = config["paths"]["proteins"]
    PRED_PATH = config["paths"]["prank"]
    POCKET_PATH = config["paths"]["pockets"]
    PAE_PATH = config["paths"]["pae"]
    OUT_PATH = config["paths"]["results"]
    PLOT_DIR = config["paths"]["plots"]
    DB_PATH = config["paths"]["db_files"]

    PLDDT_DATA = {}
    POCKET_DATA = {}
    NUM_DICT = {}

    # Create output directories if it does not exist
    if not os.path.isdir(OUT_PATH):
        os.mkdir(OUT_PATH)

    if not os.path.isdir(PLOT_DIR):
        os.mkdir(PLOT_DIR)

    # Evaluate predicted binding sites for each species
    for species in SPECIES:
        print(f"Evaluating predicted binding sites for {species}")  
        # Define paths for the current species  
        pae_path = f"{PAE_PATH}/{species}"
        struc_path = f"{STRUC_PATH}/{species}"
        pred_path = f"{PRED_PATH}/{species}"
        pocket_path = f"{POCKET_PATH}/{species}"
        out_file = f"{OUT_PATH}/bind_dict/{species}_pockets.pkl"
        tm_dict = pickle.load(open(f"{STRUC_PATH}/{species}_tm_info.pkl", "rb"))

        up_df = pd.read_csv(f"{DB_PATH}/{species}/uniprot_df.tsv.gz", delimiter="\t", index_col="Entry")
        fragments = list(up_df[~up_df["Fragment"].isna()].index)

        #Process all files in parallel and save results for plotting
        results, n_prot, plddt_all = process_files_in_parallel(struc_path, pred_path, TH, pocket_path, species, max_workers=THREADS)
        PLDDT_DATA.update(plddt_all)
        POCKET_DATA.update({pocket["pocket_id"]: pocket for pockets in results.values() for pocket in pockets})

        # Filter pockets based on pLDDT and PAE thresholds and save for further analysis
        n_pockets = 0
        for prot, pockets in results.items(): 
            pockets = [pocket for pocket in pockets if pocket["plddt"] > PLDDT_TH and pocket["pae"] < PAE_TH]
            results[prot] = pockets
            n_pockets += len(pockets)
        print(f"Found {n_pockets} pockets for {species}")
        NUM_DICT[species] = {"n_prot": n_prot, "n_pockets": n_pockets}
        with open(out_file, "wb") as f:
            pickle.dump(results, f)
    
    # Extract pLDDT and PAE data for binding sites
    pocket_df = pd.DataFrame.from_dict(POCKET_DATA, orient='index')
    pocket_df.rename(columns={"plddt": "pLDDT", "pae": "PAE"}, inplace=True)
    pocket_info = []
    plddt_list = []
    pae_list = []
    for species in SPECIES: 
        # Filter data for the current species
        tmp = pocket_df[pocket_df["species"] == species] 
        # Append pLDDT and PAE data to lists for statistical testing
        plddt_list.append(tmp["pLDDT"].to_numpy())
        pae_list.append(tmp["PAE"].to_numpy())
        
        plddt_median = tmp["pLDDT"].median()
        pae_median = tmp["PAE"].median() 
        # Obtain number of pockets before filtering, below pLDDT, and above PAE thresholds
        num_before = len(tmp)
        tmp_plddt = tmp[(tmp["pLDDT"] <= PLDDT_TH)]
        num_plddt = len(tmp_plddt)   
        tmp_pae = tmp[(tmp["PAE"] >= PAE_TH)]
        num_pae = len(tmp_pae)
        n_prot = NUM_DICT[species]["n_prot"]
        n_pockets = NUM_DICT[species]["n_pockets"]
        pocket_info.append({"species": species, "n_prot": n_prot, "n_pockets": n_pockets, "pocket/protein": round(n_pockets/n_prot, 3),
                            "plddt_med": plddt_median, "pae_med": pae_median, "num_before": num_before,
                            "num_below_plddt": num_plddt, "num_above_pae": num_pae})
    info_df = pd.DataFrame.from_records(pocket_info)

    # Statistical testing for pLDDT and PAE data
    kruskal_plddt = kruskal(*plddt_list)
    kruskal_pae = kruskal(*pae_list)
    print(f"Kruskal-Wallis test for pLDDT: {kruskal_plddt.statistic}, p-value: {kruskal_plddt.pvalue}")
    print(f"Kruskal-Wallis test for PAE: {kruskal_pae.statistic}, p-value: {kruskal_pae.pvalue}")
     
    # Plot pLDDT and PAE data for binding sites as violin plots
    dpi = 300
    fig_width = 7.5
    fig_height = 4.7
    fig, axes = plt.subplots(1, 2, figsize=(fig_width, fig_height), dpi=dpi)
    thresholds = {'pLDDT': PLDDT_TH, 'PAE': PAE_TH}
    subplot_labels = ["A", "B"]
    for ax, label, (prop, threshold) in zip(axes, subplot_labels, thresholds.items()):
        sns.violinplot(x='species', y=prop, data=pocket_df, ax=ax, inner='box', color='lightgray', linecolor='black')
        ax.axhline(y=threshold, color='firebrick', linestyle='--', linewidth=1.5, label=f'Threshold = {threshold}')
        ax.set_title(label, loc='left', fontsize=16, fontweight='bold')
        ax.set_xlabel('Species', fontsize=10)
        ax.set_ylabel(prop, fontsize=10)
        ax.tick_params("x", rotation=45)
        if prop == 'pLDDT':
            ax.legend(loc="lower right", fontsize=10)
        if prop == 'PAE':
            ax.legend(loc="upper right", fontsize=10)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    
    fig.savefig(f"{PLOT_DIR}/pocket_plddt_pae.tif", dpi=300)

    # Plot pLDDT data for full proteins 
    plddt_df = pd.DataFrame.from_dict(PLDDT_DATA, orient='index')
    
    # Define plot parameters
    fig_width = 7.5
    fig_height = 8.75
    n_row = 4
    n_col = 3
    fig, axes = plt.subplots(n_row, n_col, figsize=(fig_width, fig_height), dpi=dpi)

    prop = "prot_plddt"
    prop_name = "pLDDT"
    row = 0
    col = 0

    # Iterate through species and plot violin plots
    med_dict = {}
    for species in SPECIES:
        tmp = plddt_df[plddt_df["species"] == species]
        
        # Filter data for "With BS" and "Without BS"
        with_bs = tmp[tmp["has_pockets"] == True]
        without_bs = tmp[tmp["has_pockets"] == False]
        
        med_dict[species] = {
            "with_bs": with_bs[prop].median(),
            "without_bs": without_bs[prop].median()
        }
        
        # Combine data for plotting
        combined_data = pd.concat([
            pd.DataFrame({"Category": "With BS", prop: with_bs[prop]}),
            pd.DataFrame({"Category": "Without BS", prop: without_bs[prop]})
        ])
        sns.violinplot(
            x="Category", y=prop, data=combined_data, ax=axes[row, col],
            hue="Category", palette=["lightgreen", "orangered"],
            inner="box"
        )
        axes[row, col].set_title(species, fontsize=12, fontweight="bold")
        axes[row, col].set_xlabel("")
        axes[row, col].set_ylabel(prop_name, fontsize=10)
        
        if col == n_col - 1:
            row += 1
            col = 0
        else:
            col += 1

    axes[3, 2].set_visible(False) 
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(f"{PLOT_DIR}/prot_plddt.tif", dpi=300)

    info_df["prot_with_plddt"] = [med_dict[species]["with_bs"] for species in info_df["species"]]
    info_df["prot_without_plddt"] = [med_dict[species]["without_bs"] for species in info_df["species"]]
    info_df.to_csv(f"{PLOT_DIR}/pocket_median.csv", index=False)
