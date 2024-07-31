## DiversityInMetagenomics

This is a repository for the code scripts and synthetic and sewage data used in the master thesis project titled "Exploring Microbial Diversity in Metagenomics". Below is a description of the data objects, scripts and their functionalities.

---

## Data generation and mapping

### Python Scripts Descriptions

- **read_generator_Dirichlet1.py:**
- **read_generator_Dirichlet10.py:**
- **read_generator_sparseDirichlet.py:**
- **read_generator_sparseDirichlet2.py:**

---

### Bash Scripts Descriptions

- **RunKMA.sh:**

---

## Data processing, analytics and visualisation

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


