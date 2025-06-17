import pandas as pd 
import pickle as pkl
import yaml 
from concurrent.futures import ProcessPoolExecutor, as_completed

from Bio.PDB import PDBParser 
from Bio.PDB.SASA import ShrakeRupley


def process_pockets(prot, bind_dict, prot_dir):
    """
    Calculate the SASA for each pocket in a protein structure.
    Args:
        prot (str): Protein identifier.
        bind_dict (dict): Dictionary containing pockets for each protein.
        prot_dir (str): Directory where protein PDB files are stored.
    Returns:
        dict: A dictionary with pocket IDs as keys and their SASA and relative SASA as values.
    """
    pockets = bind_dict.get(prot)
    if not pockets: 
        return {}
    
    parser = PDBParser(QUIET=True)
    struc = parser.get_structure(prot, f"{prot_dir}/{prot}.pdb")
    
    sr = ShrakeRupley() 
    struc_sasa = sr.compute(struc, level="R")

    current_result = {}
    for pocket in pockets:
        pocket_id = pocket["pocket_id"]
        sasa = 0
        max_sasa = 0
        for res in pocket["res"]:
            sasa += struc[0]["A"][int(res)].sasa
            resname = struc[0]["A"][int(res)].get_resname()
            max_sasa += SASA_DICT.get(resname) 
        sasa = round(sasa, 3)
        rel_sasa = round(sasa/max_sasa, 3)
        current_result[pocket_id] = {"sasa": sasa, 
                                  "rel_sasa": rel_sasa}
    return current_result


def process_files_in_parallel(bind_dict, prot_dir, max_workers = None):
    result_dict = {}
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(process_pockets, prot, bind_dict, prot_dir): prot for prot in bind_dict.keys()}
        print("Start execution")
        for future in as_completed(futures):
            prot = futures[future]
            result = future.result()
            result_dict.update(result)
    return result_dict


if __name__ == "__main__":
    SASA_DICT = {"ALA": 121.0, "ARG": 265.0, "ASN": 187.0, "ASP": 187.0, "CYS": 148.0, "GLU": 214.0, "GLN": 214.0,
                "GLY": 97.0, "HIS": 216.0, "ILE": 195.0, "LEU": 191.0, "LYS": 230.0, "MET": 203.0, "PHE": 228.0,
                "PRO": 154.0, "SER": 143.0, "THR": 163.0, "TRP": 264.0, "TYR": 255.0, "VAL": 165.0} 

    # Load config file 
    with open("config.yaml", "r") as f:
        config = yaml.safe_load(f)
    N_THREADS = config["general"]["n_threads"]
    BIND_DICT_PATH = config["paths"]["bind_dict"]
    SPECIES = [species["name"] for species in config["species"]]
    STRUC_PATH = config["paths"]["proteins"]
    for species in SPECIES:
        print("Processing species:", species)
        bind_dict = pkl.load(open(f"{BIND_DICT_PATH}/{species}_pockets.pkl", "rb"))
        prot_dir = f"{STRUC_PATH}/{species}"

        results = process_files_in_parallel(bind_dict, prot_dir, N_THREADS)
        for prot, pockets in bind_dict.items():
            for pocket in pockets: 
                sasa = results.get(pocket["pocket_id"])["sasa"]
                rel_sasa = results.get(pocket["pocket_id"])["rel_sasa"]
                
                pocket["sasa"] = sasa
                pocket["rel_sasa"] = rel_sasa

        # Save updated bind dict with SASA values
        pkl.dump(bind_dict, open(f"{BIND_DICT_PATH}/{species}_pockets.pkl", "wb"))
