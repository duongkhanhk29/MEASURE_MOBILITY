## Description

This repository contains the data and code used to replicate the results from the research on intergenerational mobility, based on the Global Database on Intergenerational Mobility (GDIM 2023). The core analysis code is provided in the **R Markdown** file **Mobility_GDIM_13_Mar.RMD**, which can be used to reproduce the results of the study.

The repository includes the following files:

- **Mobility_GDIM_13_Mar.RMD**: Main R Markdown file for replication.
- **Mobility_GDIM_13_Mar - sensitive.RMD**: Script for robustness checks, to be used for alternative specifications.
- **GDIM_2023_03**: Original and processed data files for analysis.
- **GDIM_Description_2023_03**: Data description.

The country names in the dataset are shortened for data visualization. For replication, ensure all data and code are in the same folder.

### Data and Code Availability Statement

This package uses data derived from the Global Database on Intergenerational Mobility (GDIM 2023). For full details on the database, please cite the following source:

**van der Weide, Roy; Lakner, Christoph; Mahler, Daniel Gerszon; Narayan, Ambar; Ramasubbaiah, Rakesh. 2023. Intergenerational mobility around the world: A new database. Journal of Development Economics, Vol 166. [DOI](https://doi.org/10.1016/j.jdeveco.2023.103167)**

To replicate the research, download the necessary data and place it in the same directory as the provided R Markdown scripts. Ensure all required R packages, as listed in the script, are installed.

### Computational Requirements

- **Software**:
  - R (latest version recommended).
  - R Studio (for rendering the RMD files).
  - All necessary packages are specified in the script **Mobility_GDIM_13_Mar.RMD**. Install them using the `install.packages()` function.

- **Hardware**:
  - Operating System: Any modern operating system (Windows, macOS, Linux).
  - CPU: Any recent multi-core processor should be sufficient.
  - Memory: Minimum 8GB RAM recommended.
  - Disk Space: At least 1GB free space for data and results.

- **Run Time**: The analysis typically takes **less than 5 minutes** depending on your hardware configuration.

### Instructions for Data Preparation and Analysis

1. **Data Preparation**:
   - Download the required dataset from the provided source (GDIM 2023).
   - Place all data and code files in the same directory.
   - Ensure that the data file is unzipped and ready for use.

2. **Analysis**:
   - Open the **Mobility_GDIM_13_Mar.RMD** file in RStudio.
   - Follow the instructions within the script to load the data and execute the analysis.
   - For robustness checks, use the **Mobility_GDIM_13_Mar - sensitive.RMD** file following the same steps.

Both R Markdown files are designed to run without manual intervention, but ensure all necessary packages are installed before beginning the analysis.

## License

![CC-BY-NC](assets/cc-by-nc.png) This repository is under a CC-BY-NC license. Usage by commercial entities is allowed, reselling it is not.
