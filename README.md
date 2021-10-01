## LEADER_ET_AL

This repository contains matadata necessary for alignment and analysis of the human NSCLC CITEseq dataset presented in Leader, A.M., Grout, J., et al. Cancer Cell, in press (2021), also cited in Maier, B., Leader, A.M., Chen, S.T. et al. A conserved dendritic-cell regulatory program limits antitumour immunity. Nature 580, 257–262 (2020). https://doi.org/10.1038/s41586-020-2134-y

All human sequencing data is available on NCBI with BioProject ID PRJNA609924 and GEO accession GSE154826.

Alignment was performed using Cellranger v3.1.0 using feature barcoding. Feature barcode tables for the alignment are in the Leader, et al. supplemental tables.

Table S1 contains sample metadata for each sample included in the study.

Please contact andrew.leader@icahn.mssm.edu with any questions.

## Downloading the data

Data can either be downloaded automatically by running the script to reproduce the figures (see below).
Alternatively, .rd files can be downloaded using the following dropbox links:

human NSCLC scRNA & CITEseq data: https://www.dropbox.com/s/vjbide8ro5iwrfh/lung_ldm.rd?dl=1

This link will download an R data structure, the components of which contain the count matrices and cell metadata.

The data structure is called "lung_ldm" and has the following components:

1. model -> containing elements 
	1. models: a matrix with average cluster expresssion values
	2. params: a list of parameters used in the initial clustering

2. dataset -> The entire Mount Sinai 10x chromium 3' dataset presented in the paper including CITEseq data, with the following elements:
	1. umitab: raw scRNA count data of filtered cells
	2. adt_by_sample: list of raw CITEseq adt count data by sample
	3. hto_by_sample: list of raw CITEseq hto count data by sample
	4. ds: matrix of cells downsampled to 2000 UMI each
	5. cell_to_sample: array of cell to sample associations
	6. ll: log-likelihood scores for each cell mapping to each cluster
	7. ds_numis: the number of UMIs to which ds is downsampled
	8. gated_out_umitabs: raw count data for barcodes filtered during the QC filtering step
	9. counts: 3-dimensional array of samples x genes x total UMI observed per cluster
	10. samples: array of samples included in the dataset
	11. numis_before_filtering: list of arrays of number of UMIs observed per barcode in each sample prior to the filtering step
	12. max_umis: upper threshold of total UMIs per barcode used for QC filtering
	13. noise_counts: 3-dimensional array (samples x genes x cluster) of estimated # of UMI that is predicted to be attributed to noise
	14. noise_models: total average signal per sample, used as the noise component in the modified multinomial model for probabalistic classification of cells to clusters
	15. min_umis: lower threshold of total UMIs per barcode used for QC filtering
	16. avg_numis_per_sample_model: matrix of samples x clusters with values represented the average #UMI per cluster in each sample
	17. cell_to_cluster: array with cell to cluster associations
	18. alpha_noise: estimated noise fraction in each sample

## Making the figures
### Requirements

Tested on Windows 10

1. R
2. R packages (incomplete as of 10/1/21, will update shortly): 
	- gplots
	- Matrix.utils
	- mixtools
	- seriation
	- sp
	- scales
	- [scDissector](https://github.com/effiken/scDissector)

3. Downloaded and unzipped version of this repository  on a local path.

### Running the scripts in R

Assuming Leader_et_al is the local path of the repository we need to load the script files:

source("scripts/figures_main.R")`

### Output

The above referenced dropbox link will download automatically to a new /data/ directory.
Additional data files necessary to reproduce the plots will also download.

Figure will be generated in a new directory:
  - output/figures/

### Notes on specific panels
1. Figure S1A: This panel is generated during clustering, in the call to cluster() in the run_clustering.R script
2. Figures S1B, S1C: Functions to generate these plots are in the figure_s1bc.R script, but specific reproduction of these panels in the figures_main.R script has not yet been implemented.
3. Figures 7E, F and S7C, D requires downloading additional TCGA expression and mutation data. The script figure_7ef_s7cd.R performs all downloading and analysis but is not implemented inline with figures_main.R because the downloading step is time- and memory-intensive and sometimes quits unexpectedly.
4. Figures 7G-J and S7E-H analyze data from the POPLAR trial from Genentech but is not publically available.

## Clustering

### Requirements

Tested on linux LSF HPC. Due to lack of support of some of the depdendencies, the script cannot run on macOS or Windows.

1. R
2. R packages:
   - Matrix
   - Matrix.utils
   - gplots
   - seriation
   - [tglkmeans](https://github.com/tanaylab/tglkmeans)
   - [scDissector](https://github.com/effiken/scDissector)
3. Downloaded and unzipped version of this repository  on a local path.

Running the scripts in R
Assuming Leader_et_al/ is the local path of the repository, the following script will run the clustering distributedly on LSF:

source("scripts/clustering/run_clustering.r")

Note: Each run of the clustering might produce slightly different results due to different random seeds.
