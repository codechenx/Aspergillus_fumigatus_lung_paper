library(pacman)
p_load(tidyverse,phyloseq,reshape2,ape,grid,magrittr,phylosmith,parallel,DGCA,MEGENA)

metaphlan_df <- read_tsv('kraken_table.tsv')
metaphlan_df %<>% filter(grepl("Bacteria",taxon)) %>% filter(grepl("\\|s_",taxon))
metaphlan_df %<>% as.data.frame
species <- str_split(metaphlan_df$taxon, "\\|s__")
species <- sapply(species, tail, 1)
rownames(metaphlan_df) <- species
metaphlan_df$taxon <- NULL
tax_indicater <- sweep(metaphlan_df,2,colSums(metaphlan_df),`/`)
#metaphlan_df[tax_indicater<0.1]=0
meta_df <- read_csv('final_Afumigatus_censored_data.csv')
meta_df <- meta_df %>% arrange("clinic_id","cohort")
meta_df <- as.data.frame(meta_df)
meta_df <-  meta_df[c("Ifd.Nummer", "cohort","Age", "BMI","gender","clinic_id","Mutation class of allele 1","Mutation class of allele 2" )]
colnames(meta_df) <- c("Sample", "Cohort","Age", "BMI","Gender","clinic_id","Allele 1","Allele 2")
meta_df$Age <- round(meta_df$Age,1)
meta_df$Cohort<- factor(meta_df$Cohort, levels = c("neg", "before", "pos"))
rownames(meta_df) <-meta_df$Sample

taxa<- data.frame("Species" = species)
rownames(taxa) <- taxa$Species
META <- sample_data(meta_df)
TAX <- tax_table(as.matrix(taxa))
OTU <- otu_table(metaphlan_df,taxa_are_rows=TRUE)

physeq <- phyloseq(OTU, TAX, META)
physeq <- subset_samples(physeq, Cohort != "neg")
physeq <- filter_taxa(physeq,function(x) sum(x>0) >0, TRUE)
physeq <- set_sample_order(physeq, c("Cohort","clinic_id"))
SAM = sample_data(physeq) %>% data.frame(.)
stopifnot(all(SAM[SAM$Cohort =="pos",]$clinic_id == SAM[SAM$Cohort =="before",]$clinic_id))
physeq_relative_abundance <- transform_sample_counts(physeq, function(x) x / sum(x)*100 )
physeq_relative_abundance_filter <- filter_taxa(physeq_relative_abundance, function(x) sum(x > 0.1) > (0.1*length(x)), TRUE)
SAM = sample_data(physeq_relative_abundance_filter) %>% data.frame(.)
cohort <- SAM$Cohort
design_mat <- model.matrix(~0+cohort)

X<- otu_table(physeq_relative_abundance_filter)
X %<>% data.frame

set.seed(2021)
ddcor_res <-  ddcorAll(inputMat = X, design = design_mat,
  compare = c("cohortbefore", "cohortpos"),
  adjust = "perm", heatmapPlot = F,corrType = "spearman", nPerm = 999,verbose=T)

ddcor_res_significant <- ddcor_res  %>% filter(Classes != "NonSig")  %>% filter(Classes != "+/+")  %>% filter(Classes != "-/-")  %>% filter(Classes != "0/0") %>% filter(empPVals < 0.05)


set.seed(2021)
megena_res = ddMEGENA(ddcor_res_significant,
                      adjusted = FALSE,
                      nPerm=999,
                      evalCompactness = T,
                      pval_gene_thresh = 1,
                      hubPVal = 0.05,
                      modulePVal = 0.05)