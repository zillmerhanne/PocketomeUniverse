import pandas as pd 
import pickle as pkl
import os
from Bio.PDB import PDBParser 
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.ResidueDepth import ResidueDepth
from concurrent.futures import ProcessPoolExecutor, as_completed


# for pocket_id, (pocket_res, prob) in bind_dict.items():
def process_pockets(prot, bind_dict, prot_dir):
    # print(f"Get structure for {prot}")
    pockets = bind_dict.get(prot)
    if not pockets: 
        return {}
    parser = PDBParser(QUIET=True)
    struc = parser.get_structure(prot, f"{prot_dir}/{prot}.pdb")
    
    sr = ShrakeRupley() 
    # print("Calculate SASA")
    struc_sasa = sr.compute(struc, level="R")
    current_result = {}

    model = struc[0]
    # rd = ResidueDepth(model, "~/programs/msms")
    for pocket in pockets:
        pocket_id = pocket["pocket_id"]
        sasa = 0
        max_sasa = 0
        # all_depth = []
        for res in pocket["res"]:
            sasa += struc[0]["A"][int(res)].sasa
            resname = struc[0]["A"][int(res)].get_resname()
            max_sasa += SASA_DICT.get(resname) 

            # depth = rd["A",(' ', res, ' ')]
            # all_depth.append(depth[1])

        # max_depth = max(all_depth)
        # avg_depth = sum(all_depth)/len(all_depth)
        sasa = round(sasa, 3)
        rel_sasa = round(sasa/max_sasa, 3)
        current_result[pocket_id] = {"sasa": sasa, 
                                  "rel_sasa": rel_sasa}
                                  # "max_depth": max_depth,
                                  # "avg_depth": avg_depth} 
    return current_result


def process_files_in_parallel(bind_dict, prot_dir, max_workers = None):
    result_dict = {}
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(process_pockets, prot, bind_dict, prot_dir): prot for prot in bind_dict.keys()}
        print("Start execution")
        for future in as_completed(futures):
            prot = futures[future]
            result = future.result()
            # print(f"Calculated pocket SASA for {pocket_id}: {sasa}")
            result_dict.update(result)
    return result_dict


if __name__ == "__main__":
    SASA_DICT = {"ALA": 121.0, "ARG": 265.0, "ASN": 187.0, "ASP": 187.0, "CYS": 148.0, "GLU": 214.0, "GLN": 214.0,
                "GLY": 97.0, "HIS": 216.0, "ILE": 195.0, "LEU": 191.0, "LYS": 230.0, "MET": 203.0, "PHE": 228.0,
                "PRO": 154.0, "SER": 143.0, "THR": 163.0, "TRP": 264.0, "TYR": 255.0, "VAL": 165.0} 

    ORGANISMS = ["ECOLI", "YEAST", "CANAL", "ARATH", "MAIZE", "ORYSJ", "SOYBN", "DROME", "CAEEL", "MOUSE", "HUMAN"]
    for org in ORGANISMS:
        print(org)
        bind_dict = pkl.load(open(f"../../results/bind_dict/{org}_pockets.pkl", "rb"))
        prot_dir = f"../../data/proteins/{org}"
        out_dir = f"../../results/pocket_desc/{org}"
        n_workers = 150
        
        results = process_files_in_parallel(bind_dict, prot_dir, n_workers)
        for prot, pockets in bind_dict.items():
            for pocket in pockets: 
                sasa = results.get(pocket["pocket_id"])["sasa"]
                rel_sasa = results.get(pocket["pocket_id"])["rel_sasa"]
                # max_depth = results.get(pocket["pocket_id"])["max_depth"]
                # avg_depth = results.get(pocket["pocket_id"])["avg_depth"]
                pocket["sasa"] = sasa
                pocket["rel_sasa"] = rel_sasa
                # pocket["max_depth"] = max_depth 
                # pocket["avg_depth"] = avg_depth
        pkl.dump(bind_dict, open(f"../../results/bind_dict/{org}_pockets.pkl", "wb"))
