import requests
import pickle as pkl


def get_interpro_domains(uniprot_id):
    """
    Fetch domain information for a given UniProt ID from InterPro.
    """
    url = f"https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot/{uniprot_id}"
    response = requests.get(url, headers={"Accept": "application/json"})
    if response.status_code == 200:
        data = response.json()
        # print(data)
        domains = []
        for entry in data.get('results', []):
            if entry.get("metadata").get('type') == 'domain':
                for location in entry.get('proteins')[0].get("entry_protein_locations"):
                    if location["fragments"][0]["dc-status"] != "CONTINUOUS":
                        print("Found discontinuous domain")
                        print(location["fragments"])
                    domains.append({
                        "domain_id": entry["metadata"]["accession"],
                        "domain_name": entry["metadata"]["name"],
                        "start": location["fragments"][0]["start"],
                        "end": location["fragments"][0]["end"]
                    })
        return domains
    else:
        print(f"Error fetching data from InterPro: {response.status_code} for {uniprot_id}")
        return []

def map_residues_to_domains(residues, domains):
    """
    Map a list of residues to their corresponding domains.
    """
    residue_to_domain = {}
    for residue in residues:
        for domain in domains:
            if domain["start"] <= residue <= domain["end"]:
                residue_to_domain[residue] = domain["domain_name"]
                break
        else:
            residue_to_domain[residue] = "No domain"
    return residue_to_domain

# Example usage:
if __name__ == "__main__":
    organisms = ["ARATH", "ORYSJ", "MAIZE", "SOYBN", "DROME", "CAEEL", "MOUSE", "HUMAN"]
    all_spec = []
    num_bs = 0
    for org in organisms:
        print(org)
        bind_dict = pkl.load(open(f"{org}_pockets.pkl", "rb"))
        all_domains = {}
        dom_dict = {}
        all_dom = []
        # Input: UniProt ID and binding site residues
        for prot, pockets in bind_dict.items():
            uniprot_id = prot[3:]
            # print(uniprot_id)
            # Step 1: Retrieve domain information
            domains = get_interpro_domains(uniprot_id)
            all_domains[prot] = domains
            for pocket in pockets: 
                bs_res = pocket["res"]
                pocket_id = pocket["pocket_id"]
                dom_set = []
                if not domains:
                    dom_entry = {"No domain"}
                    residue_to_domain_mapping = {}
                else:
                    # Step 2: Map residues to domains
                    residue_to_domain_mapping = map_residues_to_domains(bs_res, domains)
        
                    # Step 3: Save domains 
                    dom_entry = set(list(residue_to_domain_mapping.values()))
                dom_dict[pocket_id] = {"domains": dom_entry, "res_mapping": residue_to_domain_mapping}
            all_dom += list(dom_entry)
        all_spec += all_dom
        pkl.dump(dom_dict, open(f"{org}_dom_dict.pkl", "wb"))
        all_dom_new = {}
        for key, val in all_domains.items():
            all_dom_new[key] = {x["domain_name"]: {"start": x["start"], "end": x["end"], "domain_id": x["domain_id"]} for x in val}
        all_domains = all_dom_new
        pkl.dump(all_domains, open(f"{org}_all_domains.pkl", "wb"))
