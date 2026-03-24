# Multimodal analysis of human basal ganglia astrocytes

This repository contains analysis code and release artifacts supporting a multimodal study of astrocyte diversity and organization across the adult human basal ganglia (BG).  
The work integrates single-nucleus transcriptomics and chromatin accessibility with DNA methylation, three-dimensional chromatin conformation, and spatial transcriptomics to define region- and circuit-informed astrocyte states.

The analyses presented here are part of a broader Human and Mammalian Brain Atlas (HMBA) effort within the NIH BRAIN Initiative Cell Atlas Network (BICAN), and are associated with a manuscript currently under submission.

## Overview

Astrocytes play essential roles in synaptic regulation and local homeostasis, yet their diversity and organizational principles across the human basal ganglia remain incompletely defined.  
This project focuses on establishing a multimodal, region-resolved reference of BG astrocytes and relating their molecular programs to circuit context, developmental origin, and candidate regulatory features.

Key themes include:
- Integration of transcriptomic, epigenomic, and spatial modalities
- Comparison of astrocyte states across striatal gray matter, extra-striatal gray matter, and white matter
- Identification of spatially restricted astrocyte states within conserved striatal niches
- Nomination of candidate markers and regulatory elements to guide targeted functional studies

## Data sources

Analyses in this repository draw on the following datasets:
- Human basal ganglia 10x multiome (snRNA-seq + snATAC-seq); concurrently submitted
    doi: https://doi.org/10.64898/2025.12.15.694496
- Human basal ganglia snm3C-seq; currently submitted
    doi: https://doi.org/10.64898/2026.02.12.705594
- Human and non-human primate basal ganglia spatial transcriptomic data; concurrently submitted
    doi: https://doi.org/10.1101/2025.11.22.688128
- Human whole-brain snRNA-seq reference data; public accessible
    doi: https://doi.org/10.1126/science.add7046; https://cellxgene.cziscience.com/collections/283d65eb-dd53-496d-adb7-7570c7caa443
- Mouse whole-brain single-cell and single-nucleus RNA-seq reference data; public accessible
    doi: https://doi.org/10.1038/s41586-023-06812-z; https://brain-map.org/bkp/explore/abc-atlas
- Mammalian fetal whole-brain sing-nucleus RNA-seq datasets; public accessible
    doi: https://doi.org/10.7554/eLife.109659.1; https://github.com/mtvector/scANTIPODE

Several datasets are generated as part of concurrently submitted companion manuscripts and are not yet publicly released.  
All datasets will be made publicly available in appropriate repositories upon publication.

## Repository contents

This repository includes:
- Scripts for multimodal integration, visualization, and comparative analyses

Raw sequencing data and large intermediate files are not hosted directly in this repository.

## Code availability and reproducibility

Code is provided for transparency and reproducibility of the analyses described in the associated manuscript.  
Execution of some workflows may require access to large datasets and high-performance computing resources.

## Citation

doi: https://doi.org/10.64898/2025.12.19.695583

## Contact

For questions or issues related to this repository, please contact the corresponding authors listed in the manuscript.
