# racipeCRN

[![R-CMD-check](https://github.com/76aidatill/racipeCRN/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/76aidatill/racipeCRN/actions/workflows/R-CMD-check.yaml)

*Random Circuit Perturbation for Chemical Reaction Networks*

racipECRN implements <ins>Ra</ins>ndom <ins>Ci</ins>rcuit <ins>Pe</ins>rturbation (RACIPE) for chemical reaction networks (CRNs) with mass-action kinetics. Given only the CRN's topology, it generates an ensemble of network models by randomly generating parameters. Statistical analysis of simulated results from these models reveals the basin of attraction for different states of the CRN and how they change based on extrinsic noise. The format used for the CRN topologies is Systems Biology Markup Language (SBML) or text files. See below for specifics.

## Installation ##
This is an R package. To install the development version:

```
library(devtools)
install_github("76aidatill/racipeCRN")
```

## Format for .tpo files ##
As an option besides SBML, racipeCRN allows for the creation of text files, referred to here as .tpo files, to read into the main simulation function. In these .tpo files, each line is a reaction, the reactants and products are separated by a "->" (irreversible) or "<->" (reversible) delimiter, and the stoichiometry for each reactant/product species is written as "stoichCoef speciesName" with a space separation. Like the SBML format, this system allows no products or no reactants to be specified. Below is an example .tpo text following these rules:

```
1 A -> 1 B
2 B <-> 1 C
-> 1 A
1 C ->
```
