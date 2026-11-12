# Diel_interlinked_networks

Supplement for Montoya-Bustamante et al. "Diel variation in the interlinked patterns of coupled antagonistic-mutualistic networks is associated with limited perturbation spread", Journal of Ecology.

## Instructions

If you want to reproduce the analysis and figures, first, open the `diel_interlinked.Rproj` file. Then, open any of the scripts in the "Code" folder. Scripts are independent, so there is no need to run them in any specific order. However, the way they were run in the original analysis is indicated by the number sequence in the script names.

To ensure reproducibility, this repository includes:

- A subset of the raw data is included to allow users to run and test all scripts.

- All processed data (resulting from using the full raw data) to use within the statistical models and reproduce the figures.

Raw data can be requested from the corresponding author under reasonable request.

## List of folders and files

See further info in each section.

1.  Code (folder) -\> folder containing the scripts to analyse and generate the figures of the paper. To fully reproduce figures, do not use the data files within the "Data/Processed obejcts/subset" folder.

    - `00_helping_functions.R` -\> script with a series of functions created to analyse and facilitate processing multilayer networks. Their functionality is described below:
      - `edge_list_to_bipartite()`: transforms an edgelist into a biadjacency matrix.

      - `reorder_vector()`: reorders vectors to create interlayer links in a cycle.

      - `as_weighted_biadjacency_matrix()`: transforms edge list into biadjacency matrix

      - `pref_Shannon()`: calculates weighted degree.

      - `generate_null_models()`: generates null networks based on Sauve et al. (2016)

      - `calculate_partial_cor()`: facilitates the calculation of partial correlations.
      
    - `00_diel_floquet.R` -\> Script with functions to test stability and perturbation spread using Floquet theory. Contains explanation of the functions and outputs inside the script.
    - `01_Data_preparation.R` -\> Script to create the processed objects (`.rds`). This includes all calculations of unweighted and weighted degree, correlation in consumer diversity, modularity of observed and null networks, and plant composition. To not overwrite the actual processed data used in subsequent analyses, this scripts creates new processed data within the folder "subset".
    - `02_Figure_2.R` -\> script to reproduce Figure 2.
    - `03_Figure_3.R` -\> script to reproduce Figure 3 and associated GLMMs and model outputs.
    - `04_Figure_4.R` -\> script to reproduce Figure 4 and summary outputs from modularity and similarity in module affiliation analysis.
    - `05_Data_preparation_dual.R` -\> Script to create the processed objects (`.rds`) related to dual role assignment of visitors. This includes all calculations of unweighted and weighted degree, correlation in consumer diversity, and modularity of observed and null networks. To not overwrite the actual processed data used in subsequent analyses, this scripts creates new processed data within the folder "subset".
    - `06_Figure_S1.R` -\> script to reproduce Figure S1 and associated GLMMs and model outputs.
    - `07_Figure_S2.R` -\> script to reproduce Figure S2 and summary outputs from modularity and similarity in module affiliation analysis.
    - `08_Figure_S3.R` -\> script to reproduce Figure S3.
    - `09_Figure_5_S4_S5_S6.R` -\> script to reproduce Figures 5, S4, S5, and S6, and associated tests of stability and perturbation spread.

2.  Data (folder) -\> folder containing the data used in the analyses and figures. Those files with "_dual" indicate results from analysing networks with species with dual roles.

    - `interactions_with_inttype.csv` -\> A subset of the plant-visitor interactions sampled during multiple diel periods. Columns are self-explanatory.
    - `interactions_dualint.csv` -\> A subset of the plant-visitor interactions sampled during multiple diel periods, including a classification for visitors with dual roles. Columns are self-explanatory.
    - `plant_abundances.csv` -\> a subset of the estimated plant abundances. Columns are: (1) "ID": source network, (2) "year": sampling year, (3) "site": sampling site, (4) "date": survey date, (5)"day_period" (not used) : if survey was done on a specific period, (6) "resource": plant species, (7) "abundance": floral abundance, (8) "genus": plant genus, (9) "family": plant family.
    - Processed objects (folder) -\> folder containing all processed data resulting from `01_Data_preparation.R` and `04_Figure_4.R`.
      - `degree_unweighted.rds`, `degree_unweighted_dual.rds`, and `degree_unweighted_shared_plants.rds` (data frame) -\> Correlation in consumer diversity output using unweighted degree for strict roles, dual roles, and shared plants only, respectively. Column names are self explanatory.
      - `degree_weighted.rds`, `degree_weighted_dual.rds`, and `degree_weighted_shared_plants.rds` (data frame) -\> Correlation in consumer diversity output using weighted degree for strict roles, dual roles, and shared plants only, respectively. Column names are self explanatory.
      - `Fig_2.rds` (data frame)-\> Summary data on recorded interactions used to create Figure 2. Column names are self explanatory.
      - `map_all.rds` and `map_all_dual.rds` (list) -\> Observed modularity.
      - `map_null.rds` and `map_null_dual.rds` (list) -\> Null modularity.
      - `plant_abundance_F.csv` and `polli_red_F.csv`-\> Subsets of the original raw data, removing plant species that were present but did not interact. These files are used for the stability and perturbation spread analysis. 
      - `modularity_significance_pvalues.rds` and `modularity_significance_pvalues_dual.rds` (data frame) -\> output from modularity analysis. Column names are self explanatory.
      - `plant_list.rds` and `plant_list_dual.rds` (data frame) -\> a list of the plant species in our analysis.
      - `r_af01.rds`, `r_af15.rds`, `r_mo01.rds`, `r_mo15.rds`, `r_ni01.rds`, `r_ni01.rds` -\> Results from stability and perturbation spread analysis compared to null models randomizin the mornin (mo), afternoon (af), and night (ni) networks only.
      - `rank_shift.rds` -\> data on the ranking shift of plants used to create Fig. S3.
      - `samp_intensity3.rds` and `samp_intensity3.rds` (data frame) -\> data on the number of observed interaction records per diel period. Column names are self explanatory.
      -  subset (folder) -\> a subset of the processed objects created using the subset of raw data that has been provided. Objects and explanations correspond to the ones in the "Processed objects" folder.

3.  Figures (folder) -\> contains the exported figures.

4.  Infomap.exe -\> infomap software required to run the modularity analysis using the `infomapecology` package.

## Required packages and versions

To run all analysis you will need to install multiple packages. A few of them are not currently available in CRAN.

1.  To install infomapecology see <https://infomapecology.ecomplab.com/installation.html>

2.  To install emln see <https://ecological-complexity-lab.github.io/emln_package/>.

Find a list of all packages attached to the session Info below:

```         
 progressr_0.18.0     
 furrr_0.3.1   
 future_1.67.0  
 infomapecology_2.0.2 
 vegan_2.7-2        
 permute_0.9-8        
 rlang_1.1.7      
 ggalluvial_0.12.5 
 cowplot_1.2.0        
 attempt_0.3.1        
 lmPerm_2.1.6        
 emln_1.0.2          
 hablar_0.3.2        
 DT_0.34.0           
 Matrix_1.7-4         
 ppcor_1.1           
 MASS_7.3-65         
 emmeans_2.0.0      
 DHARMa_0.4.7       
 glmmTMB_1.1.13       
 dendextend_1.19.1   
 tnet_3.0.16         
 survival_3.8-3       
 lubridate_1.9.4      
 forcats_1.0.1        
 stringr_1.6.0       
 dplyr_1.2.0        
 purrr_1.2.0         
 readr_2.1.5         
 tidyr_1.3.1          
 tibble_3.3.0      
 ggplot2_4.0.1     
 tidyverse_2.0.0     
 scales_1.4.0       
 magrittr_2.0.4      
 performance_0.15.2 
 igraph_2.2.1  
 betapart_1.6.1  
```
