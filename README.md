# racipeCRN

[![R-CMD-check](https://github.com/76aidatill/racipeCRN/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/76aidatill/racipeCRN/actions/workflows/R-CMD-check.yaml)

*Random Circuit Perturbation for Chemical Reaction Networks*

racipECRN implements <ins>Ra</ins>ndom <ins>Ci</ins>rcuit <ins>Pe</ins>rturbation (RACIPE) for chemical reaction networks (CRNs) with mass-action kinetics. Given only the topology of the CRN, it generates an ensemble of models for the circuit by randomly generating parameters. Statistical analysis of simulated results from these models reveals the basin of attraction for different states of the CRN and how they change based on extrinsic noise. The format used for the CRN topologies is Systems Biology Markup Language (SBML).

## Installation ##
This is an R package. To install the development version:

```
library(devtools)
install_github("76aidatill/racipeCRN")
```
