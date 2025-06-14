# PocketomeUniverse

## General files 
All file paths and mappings needed for the analysis are defined in the config.yaml file<br>
All scripts needed can be found in the directory "scripts"<br>
All resulting plots can be found in the directory "plots"<br>
All other results can be found in the the directory "results"<br>

## Data for binding site comparison
The following data was downloaded:
* Proteome files from UniProt for the respective species
* Predicted protein structures from [AlphaFold database](https://alphafold.ebi.ac.uk/download)
* P2Rank predictions for AlphaFold structures from [PrankWeb](https://prankweb.cz/about)
* FoldSeek cluster information was downloaded from [AlphaFold Clusters](https://afdb-cluster.steineggerlab.workers.dev/)

## Data preparation for binding site comparison 
* Binding sites predicted by P2Rank were extracted using evaluate_bs.py
* Known binding sites and information about transmembrane domains were extracted using extract_uniprot.py
* The solvent-accessible surface area for each pocket was calculated using calc_sasa.py 

## Binding site comparison using ProBis
The respective binaries for ProBiS were downloaded from their [website](http://insilab.org/probis-algorithm/). <br> 
Results were analyzed using the following scripts: 
* analyze alignments.py for extracting the alignment scores 
* cluster_probis.py for clustering the score matrices 

## Binding site comparison using DeeplyTough 
Binding site comparison using DeeplyTough was conducted using their docker image available at [docker hub](https://hub.docker.com/r/joshuameyers/deeplytough)<br> 
Results were analyzed using dim_red_all.py 

## Further analysis/visualization 
The results were further analyzed and visualized using further_analysis.ipynb
