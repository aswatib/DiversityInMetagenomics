## Exploring Microbial Diversity In Metagenomics
#### Master Thesis Project (March 1, 2024 to August 1, 2024)

This is a repository for the code scripts and synthetic and sewage data used in the master thesis project titled "Exploring Microbial Diversity in Metagenomics". Below is a description of the data objects, scripts and their functionalities.

---

## Data Generation, Mapping and Abundance Quantification

### Python Scripts Descriptions

- **read_generator_Dirichlet_alpha1.py:** generates DNA raw reads in FASTQ format with an underlying symmetric Dirichlet distribution and a concentration parameter of 10
  
- **read_generator_Dirichlet_alpha10.py:** generates DNA raw reads in FASTQ format with an underlying symmetric Dirichlet distribution and a concentration parameter of 1
  
- **read_generator_sparseDirichlet.py:** generates DNA raw reads in FASTQ format with an underlying sparse Dirichlet distribution and a concentration parameter of 0.5
  
- **read_generator_sparseDirichlet2.py:** generates DNA raw reads in FASTQ format with an underlying sparse Dirichlet distribution and a concentration parameter of 0.01
  
- **combine_bracken_report.py:** combines all individual sewage samples' Bracken-coverted Kraken2 reports and makes a single .csv file of abundance counts based on the original Kraken2 hits (excludes the bracken estimates of counts)

---

### Bash Scripts Descriptions

#### For Synthetic data

- **bacOriginal_KMA.sh** applies kma index on a set of 40 annotated complete bacterial genomes to make the BacOriginal database and then maps the synthetic raw reads to that database, outputs .res, .fsa and .mapstat files
- **bacATG_KMA.sh** maps the synthetic raw reads to the bacdb (bacATG database from CGE's KmerFinder tool), outputs .res, .fsa and .mapstat files

- **bacATG_KMA.sh:** runs KMA on synthetic reads data and maps them to the BacATG database
- **bacOriginal_KMA.sh:** runs KMA on synthetic reads data and maps them to the BacOriginal database

#### For Sewage data

- **run_KMA_Sewage.sh:** runs KMA on sewage raw reads and maps them to the bacATG database, outputs .res, .fsa and .mapstat files
- **run_Kraken2_Sewage.sh:** runs Kraken2 on sewage raw reads and maps them to the bacterial and plasmid libraries from NCBI, outputs .report files
- **run_Bracken_Sewage.sh:** runs Bracken on the individual reports obtained from Kraken2 run for the quantification of abundances from Kraken2 hits, outputs .csv with abundance counts

---

## Data Processing, Analytics and Visualisation

### R Scripts Descriptions

#### Cleaning and R packages
- **PACKAGES.R:** installs and loads all the requisite R packages and libraries used in the project, in one go
  
- **CLEANER.R:** provides functions for cleaning of species names in an abundance dataset

#### Alpha and entropy estimators
- **ALPHA_ESTIMATORS.R:** provides functions to estimate concentration parameter values for each category in an abundance dataset, applies alpha Star and alpha Star NSB from Piga et al. as well as alpha estimation using method of moments
  
- **ENTROPY_ESTIMATORS.R:** provides functions to estimate entropy using six different methods, as given by Piga et al., Nemenmann et al., Wolpert and Wolf, Hausser and Strimmer, Chao and Shen, and Shannon (1948)

#### For Synthetic data
- **SYNTHETIC_UNMAPPED.R** analytics and visualisation for unmapped raw reads for the four sets of synthetic data which was generated using python scripts above
  
- **SYNTHETIC_MAPPED.R:** analytics and visualisation for mapped counts and alignment counts using mapstat files output from the KMA run for synthetic data

- **SYNTHETIC_ENTROPIES.R:** entropy estimations for synthetic abundance datasets using mapped counts and alignment counts before and after removing variance towers

#### For Sewage data
- **SEWAGE_KMA.R:** analytics and visualisation for mapped reads using mapstat files output from the KMA run for sewage data
- **SEWAGE_KRAKEN2.R:** analytics and visualisation for mapped reads using combined bracken output csv from the Kraken2 run for sewage data

- **SEWAGE_KMA_ENTROPIES.R:** entropy estimations and other related analytics for sewage abundance datasets using mapped counts frm **KMA** hits before and after removing variance towers and also after high-noise, low-signal taxa counts (based on SNR)
- **SEWAGE_KRAKEN2_ENTROPIES.R:** entropy estimations and other related analytics for sewage abundance datasets using mapped counts from **Kraken2** hits before and after removing high-noise, low-signal taxa counts (based on SNR)

---

### RData Object Descriptions

**Abundance Tables:**
- **kma_com.RData**: Sewage abundance counts obtained from KMA hits
- **kma_after_var.RData**: Sewage abundance counts obtained from KMA hits after removing variance towers

- **kraken2_com.RData**: Sewage abundance counts obtained from Kraken2 hits

--

**Entropy Estimates:**

- **kma_entropies.RData**: entropy values for kma_com data
- **kma_after_var_entropies.RData**: entropy values for kma_com_after_var (after removing variance towers)
- **kma_after_noise_entropies.RData**: entropy values for kma_com data after removing high-noise, low-signal taxa columns
 
- **kraken2_entropies.RData**: entropy values for kraken2_com data
- **kraken2_after_noise_entropies.RData**: entropy values for kraken2_com data after removing high-noise, low-signal taxa columns

--

**Bootstrapped Entropy Estimates:**

- **kma_bootstrap_results_combinedPiga.RData**: combined bootstrap results for KMA data processed by the PIGA method
- **kraken2_bootstrap_results_combinedPiga.RData**: combined bootstrap results for Kraken2 data processed by the PIGA method

- **kma_after_var_bootstrap_resultsSingles.RData**: bootstrap results for kma_after_var data (after removing variance towers), using counts from individual samples
- **kma_after_var_bootstrap_resultsCombined.RData**: combined bootstrap results for kma_after_var data, using aggregated counts of taxa across all samples

---

## Attributions


