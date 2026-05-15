We applied the Activity-by-Contact (ABC) model (https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction) to predict cell-type-specific enhancer–promoter interactions using pseudobulk profiles of chromatin accessibility peaks and 3D contact frequency as primary inputs. 

## Setup enrivonment & download package
The original repo from broad institute:

https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction

Or:

https://github.com/DingWB/ABC

## config.yaml
See Subclass.config.yaml & ASC.config.yaml

## Prepare biosample table
Add the paths to you ATAC bigwig, peaks & hic files to biosample table (one row for each cell type / group)

## Run ABC model
Example snakemake command to run ABC model for Basal Ganglia astrocyte subGroup:
```shell
snakemake -d ~/Projects/ASC/ABC --snakefile ~/Software/ABC/workflow/run_abc.smk --configfile ASC.config.yaml -p -j 4
```