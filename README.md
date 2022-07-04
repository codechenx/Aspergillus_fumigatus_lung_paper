# Aspergillus_fumigatus_lung_paper
This repository contains the data and scripts (Python 3.7, R 4.0) that were used in the sputum microbiome by Mirhakkak et al. (2022).

### depened R packages
- pacman(0.5.1)
- tidyverse(1.3.1)
- phyloseq(1.34)
- reshape2(1.4.4)
- ape(5.6)
- grid(4.0.5)
- magrittr(2.0.2)
- phylosmith(1.0.6)
- parallel(4.0.5)
- DGCA(2.0.0)
- MEGENA(1.3.7)
- ggpubr(0.4.0)
- caret(6.0)
- circlize(0.4.13)
- ComplexHeatmap(2.6.2)
- imbalance(1.0.2.1)
- MUVR(0.0.973)
- VSURF(1.1.0)
- Boruta(7.0.0)
- mixOmics(6.8.5)
- doParallel(1.0.16)

### dependent python packages
- pandas(1.3.1)
- pycaret(2.3.2)
- imbalanced-learn(0.7.0)
- scikit-learn(0.23.2)

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