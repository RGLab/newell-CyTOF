# Newell et al. (2012) CyTOF Analysis

## Introduction

This project replicates the analysis and results from [Newell et al. (2012)](http://www.cell.com/immunity/retrieve/pii/S1074761312000040) using the CyTOF data provided.

The CyTOF data are assumed to be in the `data/` folder. Due to the large size of the files, we do not include these files in the Git repository.

## Data Munging

Before replicating the analysis, we first munge the XML workspaces and create GatingSets. The munging scripts are writtein in R and can be found in the `munge/` folder. The numbering in the filenames indicates the order in which the scripts should be executed.