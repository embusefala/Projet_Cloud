# Projet_Cloud
## TP Snakemake
The goal of this project is the initiation to a workflow language Snakemake and the discovery of its interest compared to the construction of a rudimentary workflow built in bash. 

The project is based on the construction of a workflow in Snakemake allowing to reproduce the analysis steps conducted in the ATACseq project from the quality control of the fastq.gz sequences to the identification of regions of accessibility to the DNA in 2 biological conditions considered.
Due to the impossibility to instantiate a VM with the BioPipes image (the cloud) and to install Snakemake on the Mesocenter compute cluster, I could not properly debug my snakemake script. 

The snakemake workflow corresponds to the Snakefile and the configuration files are config.yaml and env.yaml for the working environment. 
