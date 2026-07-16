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


## Repository organization

```text
BG_Astrocyte_Multimodal/
├── snRNAseq/           # snRNA-seq downstream analyses
├── snATAC/             # snATAC-seq accessibility and CREsted sequence-model analyses
│   └── CREsted/        # CREsted model training, prediction, and interpretation scripts
├── SpatialTx/          # spatial transcriptomics downstream analyses
├── snm3C-seq/          # snm3C-seq, DNA methylation, and 3D genome analyses
└── multiomic_track/    # configuration files for locus visualization and ABC modeling
    ├── *.ini           # pyGenomeTracks configuration files
    └── ABC_model/      # ABC model configuration files
```

## Processed Data sources
```text
TBD. Should be available soon!
```

## Original Raw Data sources

```text
Analyses in this repository draw on the following datasets:
- Human basal ganglia 10x multiome (snRNA-seq + snATAC-seq); concurrently submitted.
    https://doi.org/10.64898/2025.12.15.694496
- Human basal ganglia snm3C-seq; currently submitted.
    https://doi.org/10.64898/2026.02.12.705594
- Human and non-human primate basal ganglia spatial transcriptomic data; concurrently submitted.
    https://doi.org/10.1101/2025.11.22.688128
- Human whole-brain snRNA-seq reference data; public accessible.
    https://doi.org/10.1126/science.add7046; https://cellxgene.cziscience.com/collections/283d65eb-dd53-496d-adb7-7570c7caa443
- Mouse whole-brain single-cell and single-nucleus RNA-seq reference data; public accessible.
    https://doi.org/10.1038/s41586-023-06812-z; https://brain-map.org/bkp/explore/abc-atlas
- Mammalian fetal whole-brain sing-nucleus RNA-seq datasets; public accessible.
    https://doi.org/10.7554/eLife.109659.1; https://github.com/mtvector/scANTIPODE
- Genomic Tools Atlas (Allen Institute)  https://brain-map.org/bkp/experiment/genetic-tools/genetic-tools-atlas

Several datasets are generated as part of concurrently submitted companion manuscripts and are not yet publicly released.  
All datasets will be made publicly available in appropriate repositories upon publication.
```

## Related preprocessing and analysis repositories
```text
These links are provided to help readers trace the preprocessing and primary analysis steps for each reused dataset, while this repository focuses on the BG astrocyte-specific downstream analyses, revisions, and multimodal integration.

- (1) 10x Multiome / snRNA-seq and snATAC-seq data
  Primary processing and analysis workflows for the reused 10x Multiome datasets are available in the BG Consensus Taxonomy repository:  
  https://github.com/AllenInstitute/BG_Consensus_Taxonomy

  - Human snRNA-seq preprocessing and analysis:  
    https://github.com/AllenInstitute/BG_Consensus_Taxonomy/tree/main/Human/RNA

  - Human 10x snATAC-seq preprocessing and analysis:  
    https://github.com/AllenInstitute/BG_Consensus_Taxonomy/tree/main/Human/ATAC

- (2) snm3C-seq data
  Processing and analysis workflows for the reused snm3C-seq data are available here:  
  https://github.com/DingWB/BG_snm3C-seq

- (3) Spatial transcriptomics data 
  Preprocessing and analysis workflows for the reused spatial transcriptomics data are available here:  
  https://github.com/AllenInstitute/bg-spatial-manuscript-figures
```


## Multiomic Track View

https://neomorph.salk.edu/SCMDAP/  ## "Multiomic Track View" button


## Code availability and reproducibility

Code is provided for transparency and reproducibility of the analyses described in the associated manuscript.  
Execution of some workflows may require access to large datasets and high-performance computing resources.


## Citation

doi: https://doi.org/10.64898/2025.12.19.695583

## Contact

For questions or issues related to this repository, please contact the corresponding authors listed in the manuscript.
