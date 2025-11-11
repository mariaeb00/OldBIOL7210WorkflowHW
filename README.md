# OldBIOL7210WorkflowHW
This was an assignment to create a workflow using tools we had discussed in class which I copied over from my school git account
# 🧬 Workflow Assignment BIOL7210 🧬
This repo demonstrates a workflow that can complete various tasks relavent to computational genomics, specifically when looking at bacterial genomes.
There are four modules to this workflow:

1. Fastp: quality control on the raw data by filtering and trimming reads as necessary.
2. SEKSA: genome assembly with the high quality reads resulting from fastp into contigs.
3. QUAST: evaluate the quality of the previously made assemblies
4. MLST: genotyping of the previously made assemblies.

All data will be provided in the data folder within the repo. (I also uploaded them incorrectly the 1st time so they are within the main repo as well as the data folder)
- *Note to user: the data in the folder is viral genetic information, so on the MLST module you should receive an output of '-', but if you were to run this workflow on any bacterial DNA samples it will still work and provide the housekeeping genes* 

* If you would like to try bacterial samples please use the following commands or test your own :) 
```
  curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR199/000/SRR1993270/SRR1993270_*.fastq.gz
  mv SRR1993270_*.fastq.gz data/
```

This was workflow was originally developed with the following system:
- macOS Sequoia 15.1
- Architecture: x86_64
- Conda 25.1.1
  
And these are the tools that were use throughout the workflow:
- Nextflow version 24.10.5
- Fastp version 0.22.0
- SKESA version 2.3.0
- QUAST version 5.0.0
- MLST version 2.23.0


# Setting up the proper envrionment with all necessary tools on Mac ARM64 💻
```
softwareupdate --install-rosetta 
CONDA_SUBDIR=osx-64 conda create -n nf_env -c bioconda nextflow fastp skesa quast mlst
conda activate nf_env
```

# Setting up the proper envrionment with all necessary tools on Linux/Windows (I think but I am not sure, I use a mac) 💻
```
conda create -n nf_env -c bioconda -c conda-forge nextflow fastp skesa quast mlst
conda activate nf_env
```

# Usage 🎮
```
nextflow run nextflow.nf
```

