# DTANetPerturbeR

<h1 align="center">
Drug-Target binding Affinity based <br>Network Perturbation approach <br>for drug Repurposing
</h1>

<p align="center">
<a href="#introduction">Introduction</a> &nbsp;&bull;&nbsp;
<a href="#installation">Installation</a> &nbsp;&bull;&nbsp;
<a href="#usage">Usage</a>
</p>

## Introduction
**DTANetPerturbeR** is an R package designed for drug repurposing based on network perturbation algorithm. Its primary aim is to identify candidate drugs for a disease of interest. Given a set of genes associated with a disease provided by the user, the package constructs disease networks, identifies target module, and performs drug repurposing. Ultimately, it calculates repurposing scores and ranks the drugs accordingly.

### Overview of drug repurposing computational framework
![workflow](https://github.com/BarleyDavana/DTANetPerturbeR/assets/130750578/2509e5a4-a283-42cb-8096-931118b888a8)
**Drug repurposing through the PRS-based approach.** Step1: Construct a heterogeneous network comprising non-target proteins, target proteins, and drugs. Step2: Use DeepPurpose to predict binding affinities for drug-target interactions and assign weights to edges in the heterogeneous network. Step3: Employ the PRS-based algorithm to calculate the perturbation score for a specific drug acting on its target. Step4: Calculate the sum of PRS scores for each drug in the network.

### Package organization
```
DTANetPerturbeR/
├── LICENSE
├── README.md
├── R/           <- Contains R scripts for network construction and target module identification.
├── data/        <- Contains data files used by the package.
├── inst/        <- Contains Python scripts for PRS-based drug repurposing calculation.
├── man/         <- Contains .Rd help documentation.
├── DESCRIPTION  <- Package metadata.
└── NAMESPACE    <- Package namespace information.
```

## Installation

### Installation of DTANetPerturbeR
```
install.packages("devtools")
library(devtools)
install_github("BarleyDavana/DTANetPerturbeR")
```

### Installation of DeepPurpose
The DeepPurpose library is required for affinity prediction. To install it, you need to create and activate a new conda environment, install RDKit and Jupyter Notebook, install the descriptastorus dependency, and finally install DeepPurpose.
```
conda create -n DeepPurpose python=3.9
conda activate DeepPurpose
conda install -c conda-forge rdkit
conda install -c conda-forge notebook
pip install git+https://github.com/bp-kelley/descriptastorus
pip install DeepPurpose
```
For more details about DeepPurpose, please visit
https://github.com/kexinhuang12345/DeepPurpose#-a-deep-learning-library-for-compound-and-protein-modeling-dti-drug-property-ppi-ddi-protein-function-prediction

### Installation of ProDy
The ProDy package is required for PRS calculation. To install it, you can easily using the following command:
```
pip install prody
```
For more details about ProDy, please visit
https://github.com/prody/ProDy

