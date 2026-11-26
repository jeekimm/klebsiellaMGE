# klebsiellaMGE

Created on 2026-11-19

Exploring the links between phages, phage-plasmids, ARGs, and virulence factors in Klebsiella pneumoniae

Required data to run the scripts in this repo:
* phage predictions from PHASTEST, geNomad (v1.13.0), phispy (v5.0.10)
* plasmid predictions from MOB-suite (v3.1.9) and Plasmer (v23.04.20); with the Plasmer outputs mapped to PLSDB (v2024_05_31_v2) via minimap2
    
Scripts:
* compare_prophage_predictions.py : Compare PHASTEST, geNomad and PhiSpy prophage predictions, and produce consensus phages with enclosing span (i.e., inclusive method that retains the full genomic span identified by the associated predictions).
* reconcile_plasmids.py : Plasmid prediction reconciliation of MOB-suite and Plasmer results.