# Supplementary Code for "Precise Mechanochemical Scission of DNA Guided by Secondary Structures"

## Publication by: Johannes Hahmann, Arjuna Selvakumar, Boris N. Schüpp, Montgomery Labudda, Yuanxu Zhou, Gurudas Chakraborty, Frauke Gräter, and Andreas Herrmann

## GitHub maintained by: Boris N. Schüpp, current version 1.0, last updated 10/11/2025 



# How to install?
Download by using\
**git clone https://github.com/graeter-group/DNAHairpins** \
in the desired directory.

# Hardware requirements
The software runs on a standard computer with a sufficient amount of RAM. The Next Generation Sequencing code runs for the 
major dataset in less then 5 minutes.  

# Software requirements
**Next Generation Sequencing**: Software will run in a common python3 enviroment (tested with Python 3.13). \
**Required python packages**: matplotlib, seaborn, numpy, scipy \
**Molecular Dynamics**: Software can be run in a Linux terminal (only code with .sh) with GROMACS (tested on version 2025.2) and python3 installed. 
Any other code (.py) runs with in a common python3 enviroment (tested with Python 3.13). \
**Required python packages**: matplotlib, seaborn, pandas, numpy, scipy 

# What is included?

## Next Generation Sequencing Analysis Pipeline
Here, a full analysis pipeline for NGS data (.fastq) is provided. The raw data can be downloaded from [Zenodo](https://doi.org/10.5281/zenodo.17692058). NGSFilteringAlgorithm.py performs the analysis of NGS data based on a custom config file provided in the RunConfigurations directory. Additionally, all further analysis of sequencing data and the creation of plots are available in DistributionAnalysisHairpins.py, DistributionAnalysisNickedDNA.py, PlottingHairpins.py, and PlottingHairpinsSelfAssembly.py.

## Molecular Dynamics
Here, the full code to set up the simulations on an HPC environment is provided. In SimulationSetup/, the PDB files and required input files are located. In DataExtraction/, the scripts for processing trajectories are found. The raw trajectory data can be obtained from [Zenodo](https://doi.org/10.5281/zenodo.17692058). Additionally, in Analysis/, the Python code for further analysis of end-to-end distances and bond forces, as well as the code for visualization, is provided.

# Questions
For any questions regarding the provided software please contact Boris N. Schüpp (boris.schuepp@mtl.maxplanckschools.de).
