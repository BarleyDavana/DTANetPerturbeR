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
```
install.packages("devtools")
library(devtools)
install_github("BarleyDavana/DTANetPerturbeR")
```
