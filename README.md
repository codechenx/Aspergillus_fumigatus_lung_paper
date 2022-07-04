# Aspergillus_fumigatus_lung_paper
This repository contains the data and scripts (Python 3.7, R 4.0) that were used in the sputum microbiome by Mirhakkak et al. (2022).

### depened R packages
- pacman
- tidyverse
- phyloseq
- reshape2
- ape
- grid
- magrittr
- phylosmith
- parallel
- DGCA
- MEGENA
- ggpubr
- caret
- circlize
- ComplexHeatmap
- imbalance
- MUVR
- VSURF
- Boruta
- mixOmics
- doParallel

### dependent python packages
- pandas
- pycaret
- imblearn
- sklearn

### data
- growth_rate.csv --> the growth rate data of aspergillus_fumigatus for fig4B.
- min_fva_whole_v6_fastcc_95_before_kraken_final.csv -->  lower flux ranges data to support fungal growth simulated with FVA on MAMBO derived media before A. fumigatus confirmed colonization.
- max_fva_whole_v6_fastcc_95_before_kraken_final.csv --> the upper flux ranges to support fungal growth simulated with FVA on MAMBO derived media before A. fumigatus confirmed colonization.
- min_fva_whole_v6_fastcc_95_pos_kraken_final.csv --> the lower flux ranges to support fungal growth simulated with FVA on MAMBO derived media before A. fumigatus confirmed colonization.
- max_fva_whole_v6_fastcc_95_pos_kraken_final.csv --> the upper flux ranges to support fungal growth simulated with FVA on MAMBO derived media after A. fumigatus confirmed colonization.
- rxns4bp_220327.csv --> the selected reaction for fig4C.
- rxnIDsSignKoModules_220318 --> reaction id and subsystem etc data.
- kraken_table.tsv --> kraken results
- meta_data.csv --> group data

### src
- fig3_network.R --> the script to get DGCA network for drawing fig3.
- fig4B.r --> the script to create fig4B.
- fig4C.r --> the script to create fig4C.
- fig4D.r --> the script to create fig4D.
- 0.feature_selection.r 1.model_selection.py, --> a example to do feature selection as described in method.