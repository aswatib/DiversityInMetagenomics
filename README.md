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
- **bacATG_KMA.sh** maps the synthetic raw reads to the existing bacdb (bacATG database), outputs .res, .fsa and .mapstat files

- **bacATG_KMA.sh:** runs KMA on synthetic reads data and maps them to the BacATG database
- **bacOriginal_KMA.sh:** runs KMA on synthetic reads data and maps them to the BacOriginal database

#### For Sewage data

- **run_KMA_Sewage.sh:** runs KMA on sewage raw reads and maps them to the bacATG database, outputs .res, .fsa and .mapstat files
- **run_Kraken2_Sewage.sh:** runs Kraken2 on sewage raw reads and maps them to the bacterial and plasmid libraries from NCBI, outputs .report files
- **run_Bracken_Sewage.sh:** run Bracken on the individual reports obtained from Kraken2 run for the quantification of abundances from Kraken2 hits, outputs .csv with abundance counts

---

## Data Processing, Analytics and Visualisation

### R Scripts Descriptions

- **CLEAN.R:**

- **ALPHA_ESTIMATORS.R:**
- **ENTROPY_ESTIMATORS.R:**

- **SYNTHETIC_UNMAPPED.R**
- **SYNTHETIC_MAPPED.R:**

- **SEWAGE_KMA.R:**
- **SEWAGE_Kraken2.R:**

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


