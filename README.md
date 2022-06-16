# Transcriptomic and Metabolomic Response to High Light in the Charophyte Alga *Klebsormidium nitens*.

## Authors

- [@sandungxs](https://github.com/sandungxs)
- [@fran-romero-campero](https://github.com/fran-romero-campero)

## License

[GNU GPLv3](https://choosealicense.com/licenses/gpl-3.0/)

## Summary

This repository contains all the bioinformatics codes used for the realization of the TFM: **Transcriptomic and Metabolomic Response to High Light in the Charophyte Alga Klebsormidium nitens.**. They are mainly divided into two blocks according to the omics under study and the biological validation used:

- **Transcriptomic**. This folder contains all the analyses related to the study of *Klebsormidium nitens* RNA under high light (HL) stress. These include:

  - **MARACAs**. Results generated by the MARACAS software for the Klebsormidium nitens species under high light intensity for 3 hours. Outstanding among these are... 
  
  - **RNAseq**. Downstream RNAseq pipeline including global transcriptomic statistics, differential gene expression analysis and a specific gene analysis for interesting pathways.

  - **EnrichGO**. GO term enrichment for activated and repressed genes based on *Arabidopsis thaliana* orthology.

  - **Saturation_Study**. To demonstrate that using 2 replicates in this RNAseq study was sufficient, a saturation study was performed for both genome coverage and DEGs starting from different numbers of reads. In addition, a saturation study performed on an analogous *K. nitens* study performed with 3 replicates was included.

  - **Streptophyte_Unique_Genes**. To determine the set of *K.nitens* genes that are unique to Streptophyta, an orthology study was performed using various plant species and sequence similarity tools such as BLAST.





- **Metabolomic**. This folder contains all the analyses related to the study of *Klebsormidium nitens* metabolome under high light (HL) stress. These include:

  - **Downstream_analysis**. Downstream metabolomic pipeline including global metabolomic statistics, differential metabolite abundance analysis and a specific metabolite analysis for interesting pathways.
  
  - **EnrichKEGG**. KEGG pathways enrichment for activated and repressed metabolites based on *Chlorella variabilis* orthology.





- **Molecular_Biology_Validation**. This folder contains all the analyses related to the biological validation  applied to *Klebsormidium nitens* under high light (HL) stress. These include:

  - **Cell_Elongation**. Cell size measured with the software ImageJ in confocal microscopy images for each condition.
  
  - **PAM**. Chlorophyll fluorescence measured by Pulse-Amplitude-Modulation (PAM) Fluorometry for each condition.
