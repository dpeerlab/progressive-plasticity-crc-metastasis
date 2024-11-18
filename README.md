# Progressive Plasticity in Colorectal Cancer Metastasis

This repository accompanies the study *"Progressive Plasticity During Colorectal Cancer Metastasis"* published in *Nature* (available [here](https://www.nature.com/articles/s41586-024-08150-0)), aiming to reproduce key figures and analyses from the paper.

### Project Overview
The study investigates the progressive plasticity of cellular states during the metastasis of colorectal cancer. Using single-cell RNA sequencing data from primary tumors and metastases, we uncover dynamic cellular state transitions, highlighting specific lineage and state shifts associated with metastatic progression. 

### Repository Structure

- **data/**: Contains required data files:
  - **h5ads/**: AnnData files, such as `Tumor.h5ad`, `Epithelial.h5ad`, etc.
  - **tables/**: Supplementary tables in `.xlsx` format.
  - **other/**: Additional annotation and enrichment files.

- **notebooks/**: Jupyter notebooks for data download, preprocessing, and reproducing figures:
  - `download_data.ipynb`: Guide for downloading data directly from AWS S3.
  - `Figure_X.ipynb`: Notebooks for reproducing figures in the paper.
  
- **src/**: Source code modules organized by functionality, including utilities for data preprocessing (`pp`), plotting (`pl`), and label transfer or state analysis (`tl`).

### Data Access
The processed H5AD data for reproducing this analysis is hosted on AWS S3 at:
```
s3://dp-lab-data-public/progressive-plasticity-crc-metastasis
```
You can download directly using the following links:
```
https://dp-lab-data-public.s3.us-east-1.amazonaws.com/progressive-plasticity-crc-metastasis/h5ads/All.h5ad
https://dp-lab-data-public.s3.us-east-1.amazonaws.com/progressive-plasticity-crc-metastasis/h5ads/Epithelial.h5ad
https://dp-lab-data-public.s3.us-east-1.amazonaws.com/progressive-plasticity-crc-metastasis/h5ads/KG146_Organoids.h5ad
https://dp-lab-data-public.s3.us-east-1.amazonaws.com/progressive-plasticity-crc-metastasis/h5ads/KG146_Tumor.h5ad
https://dp-lab-data-public.s3.us-east-1.amazonaws.com/progressive-plasticity-crc-metastasis/h5ads/KG146_Tumor_Mapping_Reference.h5ad
https://dp-lab-data-public.s3.us-east-1.amazonaws.com/progressive-plasticity-crc-metastasis/h5ads/KG146_shPROX1_Knockdown.h5ad
https://dp-lab-data-public.s3.us-east-1.amazonaws.com/progressive-plasticity-crc-metastasis/h5ads/KG150_Tumor.h5ad
https://dp-lab-data-public.s3.us-east-1.amazonaws.com/progressive-plasticity-crc-metastasis/h5ads/KG182_Tumor.h5ad
https://dp-lab-data-public.s3.us-east-1.amazonaws.com/progressive-plasticity-crc-metastasis/h5ads/KG183_Tumor.h5ad
https://dp-lab-data-public.s3.us-east-1.amazonaws.com/progressive-plasticity-crc-metastasis/h5ads/Non-Tumor_Epithelial.h5ad
https://dp-lab-data-public.s3.us-east-1.amazonaws.com/progressive-plasticity-crc-metastasis/h5ads/Tumor.h5ad
https://dp-lab-data-public.s3.us-east-1.amazonaws.com/progressive-plasticity-crc-metastasis/h5ads/Untreated_Epithelial.h5ad
https://dp-lab-data-public.s3.us-east-1.amazonaws.com/progressive-plasticity-crc-metastasis/h5ads/Wang_etal_Tumor.h5ad
https://dp-lab-data-public.s3.us-east-1.amazonaws.com/progressive-plasticity-crc-metastasis/h5ads/Wang_etal_s1231_Tumor.h5ad
```

### Installation

To install the required dependencies, ensure you have Python 3.8 or higher and use the `pyproject.toml`:
```bash
pip install .
```

For specific package versions, review the `pyproject.toml`.

### Quickstart

1. **Download Data**: Start with `download_data.ipynb` to load data files from AWS.
2. **Run Notebooks**: Open notebooks in the `notebooks` directory to generate individual figures.

---

This README provides a basic overview of the repository's contents and is designed to support the reproducibility of key analyses and findings from the paper.
