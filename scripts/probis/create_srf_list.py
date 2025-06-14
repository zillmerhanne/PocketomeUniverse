import pandas as pd
import random as rd
import pickle as pkl
import subprocess
import os


def write_pockets(in_path, prot, pockets):
    # selection string must have the following format: "[:A and (114,116,151,153,154,181,206,22,24,254,256,259,27,280,281,282,31,55,57,58,60,89,91)]"
    for num, pocket in enumerate(pockets):
        if pocket["pae"] > 10 or pocket["plddt"] < 70 or pocket["prob"] < th:
            continue
        # pocket = sorted(pocket, key=lambda x: x[1])
        sel_string = f'"[:A and ('
        for res in pocket["res"]: 
            sel_string += str(res) + ","
        sel_string = sel_string[:-1] + ')]"'
        print(sel_string)
        # current_id = create_new_id(prot, num) 
        # subprocess.run(["gunzip", f"{in_path}/{prot}-F1-model_v4.pdb.gz"], shell=True, check=True)
        # subprocess.run(["mv", f"{in_path}/{prot}-F1-model_v4.pdb", f"{in_path}/{prot}.pdb"], shell=True, check=True)
        probis_command = [f'./probis -extract -f1 {in_path}/{prot}.pdb -c1 A -verbose -srffile {out_path}/{prot}_{num}.srf -motif {sel_string} -longnames']
        print(probis_command)
        subprocess.run(probis_command, shell=True, check=True)
        # subprocess.run(["rm", f"{in_path}/{prot}-F1-model_v4.pdb.gz"]
        # srf_files[f"{prot}_{num}"] = f"{out_path}/{current_id}.srf A\n"
        srf_files.append(f"srf/{prot}_{num}.srf A\n")
    print(f"Created {num+1} pocket surfaces for {prot}")
    return 0


def write_protein_srf(in_path, prot):
    # current_id = create_new_id_prot(prot)
    probis_command = [f'./probis -extract -f1 {in_path}/{prot}.pdb -c1 A -verbose -srffile {out_path}/{prot}.srf -longnames -out {org}']
    print(probis_command)
    subprocess.run(probis_command, shell=True, check=True)
    # srf_files[f"{prot}"] = f"{out_path}/{current_id}.srf A\n"
    srf_files.append(f"{out_path}/{prot}.srf A\n")
    print(f"Created surface for {prot}")
    return 0


if __name__ == '__main__':
    get_protein_surfaces = False
    ORGANISMS = ["ECOLI", "YEAST", "CANAL", "ARATH", "ORYSJ", "MAIZE", "SOYBN", "DROME", "CAEEL", "MOUSE", "HUMAN"]
    for org in ORGANISMS:

        prot_dir = f"../../pm_clustering/data/proteins/{org}"
        pocket_dict = pkl.load(open(f"../../pm_clustering/results/bind_dict/{org}_pockets.pkl", "rb"))
        # pae_dict = pickle.load(open(f"../../pm_clustering/results/pocket_desc/{org}/pae.pkl", "rb")) 
        # plddt_dict = pickle.load(open(f"../../pm_clustering/results/pocket_desc/{org}/plddt_dict.pkl", "rb"))
        out_path = f"{org}/srf"
        if not os.path.isdir(out_path):
            os.mkdir(out_path)
        id_file = f"srf_ids"
        db_file = f"surfaces"
        srf_files = []
    
        no_pocket = []
        has_pocket = []
        num_prot = 0
        th = 0.5

        for filename in os.listdir(prot_dir):
            if filename[-3:] == "pdb": 
                num_prot += 1
                prot_id = filename.split("/")[-1][:-4]
                if get_protein_surfaces: 
                    write_protein_srf(prot_dir, prot_id)
                else:
                    # print(prot_id)
                    pockets = pocket_dict.get(prot_id)
                    # print(pockets) 
                    if pockets:
                        write_pockets(prot_dir, prot_id, pockets)
                    else:
                        no_pocket.append(prot_id)
                        print(f"Created no pocket surfaces for {prot_id}")
    
        print(f"#proteins: {num_prot}, #proteins with no pocket: {len(no_pocket)} (= {round(len(no_pocket)/num_prot, 3)})")

        with open(f"{org}/{db_file}.txt", "w") as f: 
            for filename in srf_files: 
                f.write(filename)    
