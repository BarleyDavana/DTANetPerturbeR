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
![workflow](https://github.com/BarleyDavana/DTANetPerturbeR/assets/130750578/d6373298-da82-47c1-833a-c759ecf7f6c6)
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

### Requirements
- R 4.1.0
- Python 3.9
- pip

### Installation of DTANetPerturbeR
Use the install_github function from the devtools library to download and install the DTANetPerturbeR package.
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
https://github.com/kexinhuang12345/DeepPurpose

### Installation of ProDy
The ProDy package is required for PRS calculation. To install it, you can easily using the following command:
```
pip install prody
```
For more details about ProDy, please visit
https://github.com/prody/ProDy

## Usage

### Package architecture
DTANetPerturbeR is composed of three main steps:

* `Step1：Construct Disease Network.`
* `Step2：Identify Target Module.`
* `Step3：Perform Drug Repurposing.`

![architecture](https://github.com/BarleyDavana/DTANetPerturbeR/assets/130750578/6cebb735-339f-4af4-bfb6-361496342525)

### Construct Disease Network
The **createDiseaseNetwork.R** script is used to construct a disease protein network for a given set of disease genes.
```
├── Step1
│ ├── getCommonGenes.R
│ ├── getInitialNet.R
│ ├── runRWR.R
│ ├── getEnlargedNet.R
└── createDiseaseNetwork.R <---
```

#### Usage Example

```R
# Load disease genes from file
gene_list <- read.table("gene_list.txt", header = TRUE, sep = '\t')

# Create the disease network
createDiseaseNetwork(gene_list, enlarge_network = FALSE)
```

The `createDiseaseNetwork` function accepts two parameters:
* `disease_genes`: a data frame containing disease gene names in a single column.
* `enlarge_network`: a logical parameter indicating whether to expand the network.

It initially retrieves genes shared with other diseases and utilizes them to form an initial network. If the 'enlarge_network' parameter is set to TRUE (`default is FALSE`), the function expands the network. It ultimately returns the resulting disease network.

### Identify Target Module
The **identifyTargetModule.R** script is used to identify target module in a given disease protein network and map drugs to targets within the target module.
```
├── Step2
│ ├── getCommunities.R
│ ├── mapTargets.R
│ ├── getTargetModule.R
│ ├── mapDrugs.R
└── identifyTargetModule.R <---
```

#### Usage Example

```R
# Create a disease network using createDiseaseNetwork function
disease_net <- createDiseaseNetwork(gene_list, enlarge_network = FALSE)

# Create a disease network from a text file containing network edges
disease_net <- read.table("disease_net.txt", header = TRUE)

# Identify the target module in the network
identifyTargetModule(disease_net)
```

The `identifyTargetModule` function accepts one parameter:
* `disease_net`: a data frame containing the edge list of the disease network.

It initially detects modules within the provided disease protein network, then it maps targets to proteins within the network, identifies target module, and maps drugs to targets within the target module. It ultimately produces a data frame that represents drug-target interactions specific to the identified target module in the disease protein network.

### Perform Drug Repurposing
The **runDrugScoring.R** script is used to predict the binding affinity between drugs and targets, and calculates drug repurposing scores through PRS-based algorithm.

```
├── Step3
│ ├── predict_binding_affinity.py
│ ├── calculate_PRS.py
└── runDrugScoring.R <---
```

#### Usage Example

```R
runDrugScoring(python_exe = "python", dti_file = "Target_module_dti.txt", fasta_file = "proteins.fasta", network_file = "Target_module_edges.txt", output_folder = "output_Files")
```

The `runDrugScoring` function accepts five parameters:
* `python_exe`: The path to the Python executable.
* `dti_file`: The path to the drug-target interaction file, expected to be in TSV format with two columns: drug ID and protein UniProt ID.
* `fasta_file`: The path to the protein FASTA file.
* `network_file`: The path to the protein-protein interaction network file.
* `output_folder`: The path to the output directory where results will be saved.

It initially predicts the binding affinity between drugs and targets and calculates drug ranking scores based on these predictions by running two Python scripts. It ultimately produces a data frame with the first column representing the name of the drug, the second column listing all the targets the drug interacts with within the target module, and the third column displaying the computed repurposing score.

> **Note:** For more detailed information about each function, please refer to the function comments in the respective script.
