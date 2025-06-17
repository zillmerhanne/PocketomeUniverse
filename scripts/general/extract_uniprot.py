import pickle as pkl
import pandas as pd
import numpy as np
import os
import glob
import re
import yaml
import progressbar

from libchebipy._chebi_entity import ChebiEntity
from Bio.PDB import PDBParser

from scipy.sparse.csgraph import connected_components
from scipy.sparse import csr_matrix

import matplotlib.pyplot as plt
import seaborn as sns


def get_jaccard(list1, list2):
    """
    Calculate the Jaccard index between two lists.
    Args:
        list1 (list): First list of residues.
        list2 (list): Second list of residues.
    Returns:
        float: Jaccard index rounded to two decimal places.
    """
    # As we are mainly interested in the recovered residues of the known binding sites
    # and the fraction of residues found in a transmembrane domain, 
    # we report the Jaccard index as the fraction of residues in list1 and not the uninon of both lists.
    intersection = len(list(set(list1).intersection(list2)))
    # union = (len(list1) + len(list2)) - intersection
    return round(float(intersection)/len(list1), 2)


def get_tm_info(df, prot):
    
    if type(df.at[prot, "Transmembrane"]) == float:
        return None
    
    data = df.at[prot, "Transmembrane"]
    pattern = r'TRANSMEM\s+(\d+)(?:\.\.(\d+))?;\s+/note="([^"]+)";\s+/evidence="([^"]+)"'
    
    tm_res = []
    for match in re.finditer(pattern, data):
        start = int(match.group(1))
        end = int(match.group(2)) if match.group(2) else None
        res = list(np.arange(start, end + 1)) if end else [start]
        note = match.group(3)
        evidence = match.group(4)
        evi_start = evidence.find("ECO")
        evidence = evidence[evi_start:evi_start+11] # 11 is the length of ECO code
        tm_res.append({"res": res, "note": note, "evidence": evidence})
    return tm_res 
    


def get_bs(df, path, prot, greater_than = 4):
    """
    Extract binding sites from the UniProt database for a given protein.
    Args:
        df (pd.DataFrame): DataFrame containing UniProt data.
        path (str): Path to the PDB file of the protein.
        prot (str): UniProt ID of the protein.
        greater_than (int): Minimum number of residues in a binding site to consider it valid.
    Returns:
        list: A list of dictionaries, each containing information about a binding site.
    """
    # If the binding site data is not available for the protein, return empty lists
    if type(df.at[prot, "Binding site"]) == float:
        return []
    
    data = df.at[prot, "Binding site"]
    pocket_dict = {}
    # Regular expression to match binding site information
    pattern = r'BINDING\s+(\d+)(?:\.\.(\d+))?;\s+/ligand="([^"]+)";\s+/ligand_id="([^"]+)";(?:\s+/ligand_label="(\d+)";)?\s+/evidence="([^"]+)"'
    
    # Find all matches in the data
    for match in re.finditer(pattern, data):
        start = int(match.group(1))
        end = int(match.group(2)) if match.group(2) else None
        ligand = match.group(3)
        ligand_id = match.group(4)
        ligand_id = ligand_id.split("ChEBI:")[1] if ligand_id.startswith("ChEBI:") else ligand_id
        ligand_label = match.group(5)

        evidence = match.group(6)
        evi_start = evidence.find("ECO")
        evidence = evidence[evi_start:evi_start+11]

        if not end: 
            res = [start]
        else: 
            res = list(np.arange(start, end + 1))
        
        if (ligand, ligand_label) in pocket_dict: 
            pocket_dict[(ligand, ligand_label)]["res"].extend(res)
        else:
            pocket_dict[(ligand, ligand_label)] = {
                "res": res,
                "ligand_id": ligand_id,
                "evidence": evidence
            }
    pocket_list = []
    # Post-process the binding sites to filter and format them
    for (lig, label), pocket in pocket_dict.items():
        if len(pocket["res"]) < greater_than:
            continue
        # Check if pocket is inter-chain and split it if necessary
        split_pockets, shorter = is_interchain(path, prot, pocket["res"])
        cofactor = is_cofactor(pocket["ligand_id"])
        if not shorter:
            is_inter = False
            if len(split_pockets) > 1: 
                is_inter = True
            for res in split_pockets:
                pocket_list.append({"species": species, "ligand": pocket["ligand_id"],
                                    "lig_name": lig, "res": res, "cf": cofactor,
                                    "evidence": pocket["evidence"], "found": False,
                                    "prob": 0, "is_inter": is_inter}) 
    return pocket_list


def get_chebi(chebi, target_ids):
    """
    Cofactors:
    NADP: CHEBI:25524
    Adenosine-5'-phosphate: CHEBI:37096
    Flavine adenine dinucleutide: CHEBI:24040
    Coenzyme A: CHEBI:15346

    Compounds classes:
    Carbohydrate: CHEBI:16646
    Carbohydrate derivative: CHEBI:63299 # Glycoproteins, Amino monosaccharide (CHEBI:60926),Glycoside (CHEBI:24400), Glucoside (CHEBI:24278)
    Glycan: CHEBI:167559                   # Oligosaccharides are glycans, glucan: (CHEBI:37163) is a glycan
    Disaccharide: CHEBI:36233              # Should map to Disaccharide & derivatives
    Disaccharide derivative: CHEBI:63353   # Should map to Disaccharide & derivatives
    Monosaccharide: CHEBI:35381            # Should map to Monosaccharide & derivatives
    Monosaccharide derivative: CHEBI:63367 # Should map to Monosaccharide & derivatives

    Polypeptide: CHEBI:15841
    Oligopeptide: CHEBI:25676
    Amino acid: CHEBI:33709 # should map to AA & derivatives
    AA zwitter ion: CHEBI:35238 # should map to AA & derivatives
    AA derivatives: CHEBI:83821 # should map to AA & derivatives

    DNA: CHEBI:16991 # should map to nucleic acid
    RNA: CHEBI:33697 # should map to nucleic acid
    Nucleobase: CHEBI:18282
    Nucleoside: CHEBI:33838 # Should map to nucleosides & nucleotides
    Nucleotide: CHEBI:36976 # Should map to nucleosides & nucleotides

    Lipid: CHEBI:18059

    Inorganic ion: CHEBI:36914
    Heme: CHEBI:30413
    Hetero nuclear cluster: CHEBI:33733

    Other compound
    No known compound
    """
    chebi_entity = ChebiEntity(chebi)
    return search_relations(chebi_entity, [chebi], target_ids)


def search_relations(entity, searched, target_ids):
    """
    Recursively searches through the relations of a ChEBI entity to find a target compound class or cofactor.
    Args:
        entity (ChebiEntity): The ChEBI entity to search.
        searched (list): A list of already searched ChEBI IDs to avoid cycles.
        target_ids (list): A list of target ChEBI IDs to search for.
    Returns:
        str: The target ChEBI ID if found, otherwise None.
    """
    for rel in entity.get_outgoings():
        if rel.get_target_chebi_id() in searched or rel.get_target_chebi_id() == "CHEBI:24431":
            continue
        elif "is_" not in rel.get_type():
            continue
        elif rel.get_target_chebi_id() in target_ids:
            return rel.get_target_chebi_id()
        else:
            searched.append(rel.get_target_chebi_id())
            result = search_relations(ChebiEntity(rel.get_target_chebi_id()), searched, target_ids)
            if result:
                return result
    return None


def is_interchain(path, prot, res_list, dist_th=12, greater_than = 4):
    """
    Checks if the binding site residues are inter-chain and splits them if necessary.
    Args:
        path (str): Path to the PDB file of the protein.
        prot (str): UniProt ID of the protein.
        res_list (list): List of residue indices in the binding site.
        dist_th (float): Distance threshold to consider residues as part of the same pocket.
        greater_than (int): Minimum number of residues in a pocket to consider it valid.
    Returns:
        list: A list of lists, where each inner list contains residue indices of a pocket.
    """
    parser = PDBParser()
    structure = parser.get_structure(prot,path)
    residues = [r for r in structure.get_residues()]
    pocket_res = []
    try:
        for res in res_list:
            pocket_res.append(residues[int(res)-1])
    except:
        # Skip the pocket if there's a mismatch between the residues in the UniProt database and the AF model
        print(f"Error in {prot} with residues {res_list}. Skipping...")
        return [], True

    distances = np.empty([len(res_list), len(res_list)])
    for i in range(len(pocket_res)):
        for j in range(i, len(pocket_res)):
            res1 = pocket_res[i]
            res2 = pocket_res[j]
            distances[i, j] = res1["CA"]-res2["CA"]
            distances[j, i] = res1["CA"]-res2["CA"]

    # Inter-chain pockets in homomers, i.e., pockets that are formed by different chains, are not separated in UniProt database
    # Hence, we split them based on the all-against-all distance matrix and a threshold 
    if distances.max() > dist_th:
        dist_binary = distances
        dist_binary[dist_binary <= dist_th] = 1
        dist_binary[dist_binary != 1] = 0
        dist_binary = csr_matrix(dist_binary)
        n_components, labels = connected_components(csgraph = dist_binary, directed = False, return_labels = True)
        if n_components > 1:
            divided_pockets = [[] for i in range(n_components)]
            for i, label in enumerate(labels):
                divided_pockets[label].append(res_list[i])
            divided_pockets = [pocket for pocket in divided_pockets if len(pocket) >= greater_than]
            return divided_pockets, False
    return [res_list], False


def is_cofactor(chebi):
    """
    Checks if a given ChEBI ID corresponds to a cofactor.
    Args:
        chebi (str): The ChEBI ID to check.
    Returns:
        str: The ChEBI ID of the cofactor if found, otherwise "No cofactor" or "Other cofactor".
    """
    chebi_entity = ChebiEntity(chebi)
    for rel in chebi_entity.get_outgoings():
        if rel.get_type() == "has_role" and rel.get_target_chebi_id() == "CHEBI:23357": # CHEBI ID for biological role as cofactor
            cf = search_rel_cf(chebi_entity, [chebi])
            if cf:
                return cf
            else:
                return "Other cofactor"
    return "No cofactor"


def search_rel_cf(entity, searched):

    cf_ids = ["CHEBI:25524", "CHEBI:61296", "CHEBI:61292", "CHEBI:36981", "CHEBI:176783", "CHEBI:75769", "CHEBI:83821", "CHEBI:30413", "CHEBI:15346", "CHEBI:26191", "CHEBI:36914"]
    for rel in entity.get_outgoings():
        if rel.get_type() != "is_a" and rel.get_type() != "is_conjugate_base_of" and rel.get_type() != "is_conjugate_acid_of":
            continue
        if rel.get_target_chebi_id() in cf_ids:
            return rel.get_target_chebi_id()
        elif rel.get_target_chebi_id() in searched or rel.get_target_chebi_id() == "CHEBI:24431":
            continue
        else:
            searched.append(rel.get_target_chebi_id())
            result = search_rel_cf(ChebiEntity(rel.get_target_chebi_id()), searched)
            if result:
                 return result
    return None

if __name__ == "__main__":
    # Load configuration file and define parameters
    with open("../../config/config.yaml", "r") as f:
        config = yaml.safe_load(f)

    RESULTS_PATH = config["paths"]["results"]
    PLOT_DIR = config["paths"]["pictures"]
    DB_PATH = config["paths"]["db_files"] 
    PROT_PATH = config["paths"]["proteins"]
    BIND_DICT_PATH = config["paths"]["bind_dict"]
    KNOWN_PATH = config["paths"]["known_pockets"]
    KNOWN_FILE = f"{KNOWN_PATH}/known_pockets_up.pkl"
    if not os.path.exists(KNOWN_PATH):
        os.makedirs(KNOWN_PATH)

    SEQ_DICT = pkl.load(open(f"{PROT_PATH}/seq_len_dict.pkl", "rb"))
    SPECIES = [species["name"] for species in config["species"]]

    # Define class IDs and mappings
    # These are the ChEBI IDs for different compound classes and cofactors
    # The class_map maps ChEBI IDs to human-readable class names
    class_ids = ["CHEBI:16646", "CHEBI:63299","CHEBI:167559", 
                    "CHEBI:35381", "CHEBI:63367", "CHEBI:15841",
                    "CHEBI:25676", "CHEBI:33709", "CHEBI:35238",
                    "CHEBI:83821", "CHEBI:16991", "CHEBI:33697",
                    "CHEBI:18282", "CHEBI:33838", "CHEBI:36976",
                    "CHEBI:18059", "CHEBI:36914", "CHEBI:30413",
                    "CHEBI:33733", "CHEBI:37163", "CHEBI:24867"]
        
    class_map = {"CHEBI:16646": "Carbohydrate & derivatives", "CHEBI:63299": "Carbohydrate & derivatives", "CHEBI:167559": "Glycan", "CHEBI:37163": "Glycan",
                "CHEBI:35381": "Monosaccharide & derivatives", "CHEBI:63367": "Monosaccharide & derivatives", "CHEBI:15841": "Oligo- & polypeptide",
                "CHEBI:25676": "Oligo- & polypeptide", "CHEBI:33709": "AA & derivatives", "CHEBI:35238": "AA & derivatives",
                "CHEBI:83821": "AA & derivatives", "CHEBI:16991": "Nucleic acid", "CHEBI:33697": "Nucleic acid", "CHEBI:18282": "Nucleobase/nucleoside/nucleotide",
                "CHEBI:33838": "Nucleobase/nucleoside/nucleotide", "CHEBI:36976": "Nucleobase/nucleoside/nucleotide", "CHEBI:18059": "Lipid",
                "CHEBI:36914": "Inorganic ion", "CHEBI:24867": "Inorganic ion", "CHEBI:30413": "Heme", "CHEBI:33733": "Hetero nuclear cluster",
                "Other compound": "Other compound"}
        
    cf_map = {"CHEBI:25524": "NAD(P)", "CHEBI:61296": "ADP/ATP", "CHEBI:61292": "GDP/GTP",
            "CHEBI:36981": "FAD/FMN", "CHEBI:176783": "Other cofactor", "CHEBI:75769": "Vitamin B",
            "CHEBI:83821": "AA derivative", "CHEBI:30413": "Heme", "CHEBI:15346": "Coenzyme A", 
            "CHEBI:26191": "Other cofactor", "CHEBI:36914": "Inorganic ion", "No cofactor": "No cofactor", "Other cofactor": "Other cofactor"}

    for species in SPECIES:
        print(f"Processing {species}...")
        # Define paths for the current species
        up_path = f"{DB_PATH}/{species}"
        bind_dict_path = f"{BIND_DICT_PATH}/{species}_pockets.pkl"
        struc_path = f"../../data/proteins/{species}"

        # Load UniProt data and binding site dictionary
        up_df = pd.read_csv(os.path.join(up_path, f"uniprot_df.tsv.gz"), index_col="Entry", delimiter="\t")
        bind_dict = pkl.load(open(bind_dict_path, "rb"))
        tm_dict = pkl.load(open(f"{PROT_PATH}/{species}/{species}_tm_info.pkl", "rb")) #TODO
        tm_failed = [prot for prot, tests in tm_dict.items() if tests["num_passed"] < 6] #TODO
        print(f"Number of proteins with failed TM tests: {len(tm_failed)}")

        fragments = up_df[~up_df["Fragment"].isna()].index.tolist()
        all_prots = [os.path.basename(prot).split("-")[1] for prot in glob.glob(f"{struc_path}/*pdb")]
        all_prots = [prot for prot in all_prots if prot not in fragments and SEQ_DICT.get(prot, 0) > 100 and prot not in tm_failed]

        tm_dom_dict = {}
        known_bs = {}
        lig_map = {}

        with progressbar.ProgressBar(max_value=len(all_prots)) as bar:
            for idx, prot in enumerate(all_prots):
                if prot not in up_df.index:
                    ec = "No EC"
                    known_bs[prot] = []
                    tm_dom_dict[prot] = []
                    continue

                # Get EC number for the protein
                ec = up_df.at[prot, "EC number"]
                if type(ec) == float:
                    ec = "No EC"
                elif len(ec.split(";")) > 1:
                    ec = "More than 1 EC"
                else:
                    ec = f"EC {int(ec[0])}.-.-.-"

                # Extract binding sites for the protein
                pockets = get_bs(up_df, os.path.join(struc_path, f"AF-{prot}.pdb"), prot)
                # Assign ligand classes and names to the pockets
                for pocket in pockets:
                    pocket["EC"] = ec 
                    if pocket["ligand"] not in lig_map.keys():
                        compound_class = get_chebi(pocket["ligand"] , class_ids)
                        if not compound_class or compound_class == "CHEBI:24431":
                            compound_class = "Other compound"
                        class_name = class_map.get(compound_class, "Other compound")
                    else: 
                        compound_class = lig_map.get(pocket["ligand"])
                        class_name = class_map.get(compound_class, "Other compound")
                    pocket["lig_class"] = compound_class
                    pocket["class_name"] = class_name
               
                known_bs[prot] = pockets
                tm_domains = get_tm_info(up_df, prot)
                tm_dom_dict[prot] = tm_domains
                bar.update(idx)

        for prot, pockets in bind_dict.items():
            tm_domains = tm_dom_dict.get(prot[3:], None)
            for pocket in pockets:

                prob = pocket["prob"]
                known_pockets = known_bs.get(prot[3:])
                if not known_pockets: 
                    known_pockets = []

                # If the pocket is found in the known pockets, we assign it:
                # - a ligand
                # - a class name
                # - a cofactor class
                # - a Jaccard index (overlap with known pocket)
                # If the pocket is not found in the known pockets, we assign it:
                # - "No known ligand" as compound class 
                # - "No cofactor" as cofactor class
                label = "No known ligand"
                cf_label = "No cofactor"
                chebi = None
                jacc = None
                for idx, known in enumerate(known_pockets):
                    if list(set(known["res"]) & set(pocket["res"])):
                        if chebi:
                            chebi.append(known["ligand"])
                            label = "Several ligands"
                            jacc.append(get_jaccard(known["res"], pocket["res"]))
                        else:
                            jacc = [get_jaccard(known["res"], pocket["res"])]
                            chebi = [known["ligand"]]
                            label = known["class_name"]

                        cf_label = cf_map.get(known["cf"])
                        
                        if known["found"] == False:
                            known["found"] = True
                            known["prob"] = prob 
                            known["jacc"] = get_jaccard(known["res"], pocket["res"])

                        # A known pocket can be intersect with several predicted pockets
                        # In these cases, we keep the one with the highest probability
                        elif known["found"] == True and known["prob"] < prob:
                            known["prob"] = prob
                            known["jacc"] = get_jaccard(known["res"], pocket["res"])

                if chebi and len(chebi) == 1: 
                    chebi = chebi[0]
                    jacc = jacc[0]
                # If the pocket is not found in the known pockets, we assign it a new ligand
                pocket["ligand"] = chebi
                pocket["lig_class"] = label
                pocket["jacc"] = jacc
                pocket["cf_class"] = cf_label

                # Check if pocket is found in a transmembrane domain
                # If so, add the information to the pocket
                is_tm = False 
                tm_note = [] 
                tm_evidence = [] 
                jacc = None
                if tm_domains:
                    for tm in tm_domains:
                        if list(set(tm["res"]) & set(pocket["res"])):
                            is_tm = True
                            tm_note.append(tm["note"])
                            tm_evidence.append(tm["evidence"])
                            jacc = get_jaccard(pocket["res"], tm["res"])
                pocket["is_tm"] = is_tm
                pocket["tm_note"] = set(tm_note)
                pocket["tm_evidence"] = set(tm_evidence)
                pocket["tm_jacc"] = jacc

        # Save the updated binding site dictionary and known pockets
        if not os.path.isfile(KNOWN_FILE):
            pkl.dump(known_bs, open(KNOWN_FILE, "wb"))
        else:
            known_up = pkl.load(open(KNOWN_FILE, "rb"))
            known_up.update(known_bs)
            pkl.dump(known_up, open(KNOWN_FILE, "wb"))

        pkl.dump(bind_dict, open(bind_dict_path, "wb"))