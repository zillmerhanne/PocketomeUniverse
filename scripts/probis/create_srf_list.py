import yaml
import pickle as pkl
import subprocess
import os
import glob


def write_pockets(in_path, prot, pockets):
    """
    Write the pockets of a protein to SRF files using ProBiS.
    Args:
        in_path (str): Path to the directory containing protein PDB files.
        prot (str): Protein identifier.
        pockets (list): List of pockets, each containing a list of residue indices.
    Returns:
        int: 0 if successful
    """
    # selection string must have the following format: "[:A and (114,116,151,153,154,181,206,22,24,254,256,259,27,280,281,282,31,55,57,58,60,89,91)]"
    for num, pocket in enumerate(pockets):
        sel_string = f'"[:A and ('
        for res in pocket["res"]: 
            sel_string += str(res) + ","
        sel_string = sel_string[:-1] + ')]"'
        probis_command = [f'./probis -extract -f1 {in_path}/{prot}.pdb -c1 A -verbose -srffile {out_path}/{prot}_{num}.srf -motif {sel_string} -longnames']
        subprocess.run(probis_command, shell=True, check=True)
        srf_files.append(f"{out_path}/{prot}_{num}.srf A\n")
    return 0


def write_protein_srf(in_path, prot):
    """
    Write the surface of a protein to an SRF file using ProBiS.
    Args:
        in_path (str): Path to the directory containing protein PDB files.
        prot (str): Protein identifier.
    Returns:
        int: 0 if successful"""
    probis_command = [f'./probis -extract -f1 {in_path}/{prot}.pdb -c1 A -verbose -srffile {out_path}/{prot}.srf -longnames -out {species}']
    subprocess.run(probis_command, shell=True, check=True)
    srf_files.append(f"{out_path}/{prot}.srf A\n")
    return 0


if __name__ == '__main__':
    get_protein_surfaces = False
    # Load conifgurations
    with open("config.yaml", "r") as f:
        config = yaml.safe_load(f)
    SPECIES = [species["id"] for species in config["species"]]
    STRUC_PATH = config["paths"]["proteins"]
    BIND_DICT_PATH = config["paths"]["bind_dict"]
    PROBIS_PATH = config["paths"]["probis_results"]
    for species in SPECIES:
        print(f"Processing species: {species}")
        # Set species-specific paths
        prot_dir = f"{STRUC_PATH}/{species}"
        pocket_dict = pkl.load(open(f"{BIND_DICT_PATH}/{species}_pockets.pkl", "rb"))
        out_path = f"{PROBIS_PATH}/{species}/srf"
        if not os.path.isdir(out_path):
            os.mkdir(out_path)

        db_file = f"surfaces"
        srf_files = []

        for filename in glob.glob(f"{prot_dir}/*.pdb"):
            filename = os.path.basename(filename) 
            prot_id = "-".join(filename.split("-")[:2])  # Extract protein ID from filename
            if get_protein_surfaces: 
                write_protein_srf(prot_dir, prot_id)
            else:
                pockets = pocket_dict.get(prot_id)
                if pockets:
                    write_pockets(prot_dir, prot_id, pockets)
    
        # Write surface database for probis comparison 
        with open(f"{PROBIS_PATH}/{species}/{db_file}.txt", "w") as f: 
            for filename in srf_files: 
                f.write(filename)    
