import pickle as pkl
import pandas as pd
import numpy as np
import os
import glob
import re
from libchebipy._chebi_entity import ChebiEntity
from Bio.PDB import PDBParser
from scipy.sparse.csgraph import connected_components
from scipy.sparse import csr_matrix
import progressbar
import matplotlib.pyplot as plt
import seaborn as sns


def get_jaccard(list1, list2):
    intersection = len(list(set(list1).intersection(list2)))
    union = (len(list1) + len(list2)) - intersection
    return round(float(intersection)/len(list1), 2)


def get_tm_info(df, prot):
    # TO DO
    # Include ligand label --> discriminates between different ligand entities
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
    


def get_bs(df, path, prot, greater_than = 4, shorter_than=32):
    if type(df.at[prot, "Binding site"]) == float:
        return [], {}, []
    
    data = df.at[prot, "Binding site"]
    pocket_dict = {}
    pattern = r'BINDING\s+(\d+)(?:\.\.(\d+))?;\s+/ligand="([^"]+)";\s+/ligand_id="([^"]+)";(?:\s+/ligand_label="(\d+)";)?\s+/evidence="([^"]+)"'
    
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
    ligands = {}
    lig_list = []
    for (lig, label), pocket in pocket_dict.items():
        num_res = len(pocket["res"])
        if pocket["ligand_id"] == "CHEBI:30413" or "heme" in lig.lower():
            print(f"Found heme binding site for : {prot}, with {num_res} residues")
        # Check if pocket is inter-chain 
        split_pockets, shorter = is_interchain(path, prot, pocket["res"])
        cofactor = is_cofactor(pocket["ligand_id"])
        if not shorter:
            is_inter = False
            if len(split_pockets) > 1: 
                is_inter = True
            for res in split_pockets:
                if len(res) >= greater_than and len(res) <= shorter_than:
                    pocket_list.append({"org": org, "ligand": pocket["ligand_id"], "lig_name": lig, "res": res, "cf": cofactor, "evidence": pocket["evidence"], "found": False, "prob": 0, "is_inter": is_inter})      # Check if AF model has less AAs than where the pocket is located, e.g., O60673, and that pocket is longer than a certain threshold
                    ligands[pocket["ligand_id"]] = 1
                    lig_list.append(pocket["ligand_id"])
    return pocket_list, ligands, lig_list


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

    Other compounds:
    Vitamin B: CHEBI:75769
    Vitamin E: CHEBI:33234
    Vitamin K: CHEBI:28384
    Quinone: CHEBI:36141
    Sulfide: CHEBI:26822
    Chlorophyll: CHEBI:28966
    Heteroarene: CHEBI:33833  # semaxanib (azole)
    Ester: CHEBI:35701
    Amide: CHEBI:32988
    Aldehyde: CHEBI:17478
    Alcohol: CHEBI:30879
    Thiazide: CHEBI:50264
    Organohalogen compound: CHEBI:17792
    Elemental molecular entity: CHEBI:33259 #Oxygen, carbon, etc
    Iron chelate: CHEBI:5975
    Elemental molecular entity: CHEBI:33259 #Oxygen, carbon, etc
    """
    chebi_entity = ChebiEntity(chebi)
    return search_relations(chebi_entity, [chebi], target_ids)


def search_relations(entity, searched, target_ids):
    for rel in entity.get_outgoings():
        if rel.get_target_chebi_id() in searched or rel.get_target_chebi_id() == "CHEBI:24431":
            continue
        # elif rel.get_type() != "is_a" and rel.get_type() != "is_conjugate_base_of" and rel.get_type() != "is_conjugate_acid_of" and rel.:
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
    parser = PDBParser()
    structure = parser.get_structure(prot,path)
    residues = [r for r in structure.get_residues()]
    pocket_res = []
    try:
        for res in res_list:
            pocket_res.append(residues[int(res)-1])
    except:
        print(prot, res)
        return [], True

    distances = np.empty([len(res_list), len(res_list)])
    for i in range(len(pocket_res)):
        for j in range(i, len(pocket_res)):
            res1 = pocket_res[i]
            res2 = pocket_res[j]
            distances[i, j] = res1["CA"]-res2["CA"]
            distances[j, i] = res1["CA"]-res2["CA"]

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

# Define in separate config file 
# No absoulute path should be defined within the code
RESULTS_PATH = "../../results"
PLOT_DIR = "../../results/pictures/REVISION"
DB_PATH = "../../../Pocketeome/db_files"
SEQ_PATH = "../../results/proteins/seq_len_dict.pkl"
SEQ_DICT = pkl.load(open(SEQ_PATH, "rb"))
ORGANISMS = ["ECOLI", "YEAST", "CANAL", "ARATH", "ORYSJ","MAIZE", "SOYBN", "DROME", "CAEEL", "MOUSE", "HUMAN"]
# ORGANISMS = ["ECOLI"]
ALL_POCKETS = []
KNOWN_FILE = "../../data/known_pockets/known_pockets_up.pkl"
RESULTS = []

for org in ORGANISMS:
    if not os.path.isdir(f"../../results/pictures/{org}"):
        os.mkdir(f"../../results/pictures/{org}")

    up_path = f"~/Pocketeome/db_files/{org}"
    bind_dict_path = f"../../results/all_pockets/{org}_pockets.pkl"
    out_path = f"../../results/pocket_desc/{org}"
    struc_path = f"../../data/proteins/{org}"

    up_df = pd.read_csv(os.path.join(up_path, f"{org}_UP.csv"), index_col="Entry")
    bind_dict = pkl.load(open(bind_dict_path, "rb"))
    tm_dict = pkl.load(open(f"{RESULTS_PATH}/pocket_desc/{org}/{org}_tm_info.pkl", "rb"))
    tm_failed = [prot for prot, tests in tm_dict.items() if tests["num_passed"] < 6]
    print(f"Number of proteins with failed TM tests: {len(tm_failed)}")

    if os.path.isfile(os.path.join(up_path, f"frag_df.tsv.gz")):
        fragments = list(pd.read_csv(os.path.join(up_path, f"frag_df.tsv.gz"), index_col="Entry", delimiter="\t", compression="infer").index)
    else: 
        fragments = []
    all_prots = [os.path.basename(prot)[3:-4] for prot in glob.glob(f"{struc_path}/*pdb")]
    all_prots = [prot for prot in all_prots if prot not in fragments and SEQ_DICT.get(prot, 0) > 100 and prot not in tm_failed]

    # Define class IDs and class map in config file
    # These IDs are from ChEBI and should be updated if new classes are added or changed 
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

    class_color = {"No known ligand": "#C7C7C7", "Other compound": "gray", "Several ligands": "#343534", "Monosaccharide & derivatives": "#AA0DFE", 
                "Glycan": "#800080", "Carbohydrate & derivatives": "#FF00FF", "AA & derivatives": "#00FF00", "Oligo- & polypeptides": "#109618",
                "Nucleobase/nucleoside/nucleotide": "#0000FF", "Nucleic acid": "#40E0D0", "Lipid": "#FF0000", "Inorganic ion": "#FF9DA6",
                "Hetero nuclear cluster": "#FF7F00", "Heme": "#FFFF00"}

    plant_ids = ["CHEBI:28966", "CHEBI:23044", "CHEBI:26125", "CHEBI:24279", "CHEBI:72544", 
                "CHEBI:26848", "CHEBI:6457", "CHEBI:26873", "CHEBI:22315", "CHEBI:26605",
                "CHEBI:22676", "CHEBI:23530", "CHEBI:24250"]
    
    cf_map = {"CHEBI:25524": "NAD(P)", "CHEBI:61296": "ADP/ATP", "CHEBI:61292": "GDP/GTP",
            "CHEBI:36981": "FAD/FMN", "CHEBI:176783": "Other cofactor", "CHEBI:75769": "Vitamin B",
            "CHEBI:83821": "AA derivative", "CHEBI:30413": "Heme", "CHEBI:15346": "Coenzyme A", 
            "CHEBI:26191": "Other cofactor", "CHEBI:36914": "Inorganic ion", "No cofactor": "No cofactor", "Other cofactor": "Other cofactor"}

    cf_colors = {"ADP/ATP": "#1980ed", "Coenzyme A": "#2ca02c", "FAD/FMN": "#d62728", "GDP/GTP": "#9467bd",
                "Heme": "#f0e442", "NAD(P)": "#ff9896", "Inorganic ion": "#7d3c98", "AA derivative": "#fc427b",
                "Vitamin B": "#008080", "Other cofactor": "gray", "No cofactor": "#bebdbd"}


    ec_dict = {}
    cofactor_dict = {}
    tm_dom_dict = {}
    known_bs = {}
    known_bs_found = {}
    all_ligands = {}
    all_lig_list = []
    # n_inter = 0
    print(org)

    ligands = {target:0 for target in class_map.keys()}
    ligands["Other compound"] = 0
    lig_map = {}

    with progressbar.ProgressBar(max_value=len(all_prots)) as bar:
        for idx, prot in enumerate(all_prots):
            if prot not in up_df.index:
                ec = "No EC"
                ec_dict[prot] = "No EC"
                known_bs[prot] = []
                tm_dom_dict[prot] = []
                continue

            ec = up_df.at[prot, "EC number"]

            if type(ec) == float:
                ec_dict[prot] = "No EC"
                ec = "No EC"
            elif len(ec.split(";")) > 1:
                ec_dict[prot] = "More than 1 EC"
                ec = "More than 1 EC"
            else:
                ec_dict[prot] = f"EC {int(ec[0])}.-.-.-"
                ec = f"EC {int(ec[0])}.-.-.-"
            
            pockets, ligands, ligand_list = get_bs(up_df, os.path.join(struc_path, f"AF-{prot}.pdb"), prot)
            tm_domains = get_tm_info(up_df, prot)

            for pocket in pockets:
                pocket["EC"] = ec 
                if pocket["cf"] not in cofactor_dict.keys():
                    cofactor_dict[pocket["cf"]] = 1
                elif pocket["cf"] in cofactor_dict.keys():
                    cofactor_dict[pocket["cf"]] += 1

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

            # n_inter += current_inter
            tm_dom_dict[prot] = tm_domains
            known_bs[prot] = pockets
            known_bs_found[prot] = [False] * len(pockets)
            all_ligands.update(ligands)
            all_lig_list += ligand_list
            bar.update(idx)

    ec_labels = {}
    lig_labels = {}
    cf_labels = {}
    PRED_POCKETS =  []

    for prot, pockets in bind_dict.items():
        tm_domains = tm_dom_dict.get(prot[3:], None)
        for pocket in pockets:
            PRED_POCKETS.append(pocket)

            prob = pocket["prob"]
            known_pockets = known_bs.get(prot[3:])
            if not known_pockets: 
                known_pockets = []

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
                    elif known["found"] == True and known["prob"] < prob:
                            known["prob"] = prob
                            known["jacc"] = get_jaccard(known["res"], pocket["res"])

            if chebi and len(chebi) == 1: 
                chebi = chebi[0]
                jacc = jacc[0]

            pocket["ligand"] = chebi
            pocket["lig_class"] = label
            pocket["jacc"] = jacc
            pocket["cf_class"] = cf_label

            is_tm = False 
            tm_note = [] 
            tm_evidence = [] 
            if tm_domains:
                for tm in tm_domains:
                    if list(set(tm["res"]) & set(pocket["res"])):
                        is_tm = True
                        tm_note.append(tm["note"])
                        tm_evidence.append(tm["evidence"])
            pocket["is_tm"] = is_tm
            pocket["tm_note"] = set(tm_note)
            pocket["tm_evidence"] = set(tm_evidence)

    if not os.path.isfile(KNOWN_FILE):
        pkl.dump(known_bs, open(KNOWN_FILE, "wb"))
    else:
        known_up = pkl.load(open(KNOWN_FILE, "rb"))
        known_up.update(known_bs)
        pkl.dump(known_up, open(KNOWN_FILE, "wb"))

    pkl.dump(bind_dict, open(bind_dict_path, "wb"))