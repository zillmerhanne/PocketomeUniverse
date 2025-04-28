# PocketomeUniverse

## Data preparation for binding site comparison
For prediction and processing of potential binding sites a Snakemake pipeline was implemented. 
This pipeline includes the following steps: 
1. Download of protein structures
2. Pre-processing of protein structures 
3. Pocket prediction
4. Pocket evaluation
5. Pocket feature extraction  

## Binding site comparison using ProBis
The respective binaries for ProBiS were downloaded from their [website](http://insilab.org/probis-algorithm/). 
The following commands were used to compare binding site using ProBis: 
Analysis scripts for binding site comparison via ProBis can be found in the respective sub-folder.

## Binding site comparison using DeeplyTough 
Binding site comparison using DeeplyTough was conducted using [their docker image available at docker hub](https://hub.docker.com/r/joshuameyers/deeplytough) and the following commands: 
Analysis scripts for binding site comparison via ProBis can be found in the respective sub-folder.
