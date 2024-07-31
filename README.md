## DiversityInMetagenomics

This is a repository for the code scripts and synthetic data used in the master thesis project titled "Exploring Microbial Diversity in Metagenomics". Below is a description of the data objects, scripts and their functionalities.

---

### RData Object Descriptions

**Abundance Tables:**
- **kma_com.RData**: Sewage abundance counts obtained from KMA hits
- **kma_after_var.RData**: Sewage abundance counts obtained from KMA hits after removing variance towers

- **kraken2_com.RData**: Sewage abundance counts obtained from Kraken2 hits

**Entropy Estimates:**

- **kma_entropies.RData**: entropy values for raw KMA data
- **kma_after_var_entropies.RData**: entropy values for KMA data after removing variance towers
- **kma_after_noise_entropies.RData**: entropy values for KMA data after high-noise, low-signal taxa removal
 
- **kraken2_entropies.RData**: entropy values for raw Kraken2 data
- **kraken2_after_noise_entropies.RData**: entropy values for Kraken2 data after high-noise, low-signal taxa removal

**Bootstrapped Entropy Estimates:**

- **kma_bootstrap_results_combinedPiga.RData**: combined bootstrap results for KMA data processed by the PIGA method
- **kraken2_bootstrap_results_combinedPiga.RData**: combined bootstrap results for Kraken2 data processed by the PIGA method

- **kma_after_var_data_bootstrap_resultsSingles.RData**: bootstrap results for KMA data after removing variance towers, focusing on individual samples
- **kma_after_var_data_bootstrap_resultsCombined.RData**: combined bootstrap results for KMA data after removing variance towers, focusing on aggregated counts of taxa across samples

---


