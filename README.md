# emergeR
## en Masse Evaluation of RNA Guides: R suite
In development.  
Developer: Anthony Surkov (Beal Lab, UC Davis)  
Objective: end-to-end pipeline for analysis of high-throughput RNA editing data, including processing, statistics, visualization, motif identification, and sequence candidate selection.  
More information on EMERGe: https://pubs.acs.org/doi/full/10.1021/acschembio.3c00107?casa_token=7Cc8Heg_HzwAAAAA%3Aanlj8gTCUpNP27y921NhgH0M_j4lUfDGJbS9SyyIs0cK3rap-Vn_ILNAjljLUNTOI8X_zPLonpE-Lqk  

Status:
- Main pipeline (R/) is an adaptation of disparate tools used in current analysis (https://doi.org/10.1021/scimeetings.5c11372). Currently undergoing testing.
- Additional tools (dev/) are being built to resolve data processing issues (fixed Z base not being read as C, regex preprocessing helping with homopolymeric repeats in NGS data) and generate reports for EMERGe screens.
- Planned work includes ViennaRNA utilization to predict viability of guide-target binding in an intermolecular reaction context.
