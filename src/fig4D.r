library(tidyverse)
library(circlize)
library(ComplexHeatmap)
library(ggpubr)

#Loading data
selected_rxn_df <- read_csv("../data/rxns4bp_220327.csv")

filter_df <- read_tsv("../data/rxnIDsSignKoModules_220318.tsv")
pos_max_dataset <- read.csv('../data/max_fva_whole_v6_fastcc_95_pos_kraken_final.csv')
rownames(pos_max_dataset) <- paste0("max_",pos_max_dataset$Metabolite)
pos_max_dataset$Metabolite<- NULL
pos_min_dataset <- read.csv('../data/min_fva_whole_v6_fastcc_95_pos_kraken_final.csv')
rownames(pos_min_dataset) <-paste0("min_",pos_min_dataset$Metabolite)
pos_min_dataset$Metabolite<- NULL
before_max_dataset <- read.csv('../data/max_fva_whole_v6_fastcc_95_before_kraken_final.csv')
rownames(before_max_dataset) <-paste0("max_",before_max_dataset$Metabolite)
before_max_dataset$Metabolite<- NULL
before_min_dataset <- read.csv('../data/min_fva_whole_v6_fastcc_95_before_kraken_final.csv')
rownames(before_min_dataset) <- paste0("min_",before_min_dataset$Metabolite)
before_min_dataset$Metabolite<- NULL

max_dataset <- cbind(before_max_dataset,pos_max_dataset) %>% t %>% as.data.frame()
min_dataset <- cbind(before_min_dataset,pos_min_dataset) %>% t %>% as.data.frame()



# calculate p value
max_dataset_p <- data.frame()
for (metabolite in colnames(max_dataset)){
  wilcox_df <- max_dataset[metabolite]
  wilcox_df[wilcox_df==0] <- 0.00001
  wilcox_v <- wilcox_df[[metabolite]]
  if(var(na.omit(wilcox_v)) == 0){
    max_dataset_p<- rbind(max_dataset_p, data.frame("metabolite"=metabolite, "p"=1,fold_change=1))
    next
  }
  wilcox_x <- wilcox_v[1:49]
  wilcox_y <-wilcox_v[50:98]
  fold_change <- median(na.omit(wilcox_y/wilcox_x))
  wilcoxmodel <- wilcox.test(wilcox_x,wilcox_y,paired=TRUE)
  p_value <- wilcoxmodel$p.value
  max_dataset_p<- rbind(max_dataset_p, data.frame("metabolite"=metabolite, "p"=p_value,fold_change=fold_change))
}

max_dataset_p$adjp <- p.adjust(max_dataset_p$p, "fdr")
max_dataset_p$difffc <- abs(max_dataset_p$fold_change-1)

min_dataset_p <- data.frame()
for (metabolite in colnames(min_dataset)){
  wilcox_df <- min_dataset[metabolite]
  wilcox_df[wilcox_df==0] <- 1
  wilcox_v <- wilcox_df[[metabolite]]
  if(var(na.omit(wilcox_v)) == 0){
     min_dataset_p<- rbind(min_dataset_p, data.frame("metabolite"=metabolite, "p"=1,fold_change=1))
     next
  }
  wilcox_x <- wilcox_v[1:49]
  wilcox_y <-wilcox_v[50:98]
  fold_change <- median(na.omit(wilcox_y/wilcox_x))
  wilcoxmodel <- wilcox.test(wilcox_x,wilcox_y,paired=TRUE)
  p_value <- wilcoxmodel$p.value
  min_dataset_p<- rbind(min_dataset_p, data.frame("metabolite"=metabolite, "p"=p_value,fold_change=fold_change))
}

min_dataset_p$adjp <- p.adjust(min_dataset_p$p, "fdr")
min_dataset_p$difffc <- abs(min_dataset_p$fold_change-1)

selected_rxn_final <- c("R1965","R1838","R153")
selected_rxn__final_name <- c("EC 4.1.3.27","EC 4.2.1.51","EC 1.8.1.9")

# plot boxplot
plots = NULL
for(i in 1:2){
for(scale in 1:length(selected_rxn_final)) {
        rxn_df <- filter_df %>% filter(ID_new == selected_rxn_final[scale]) 
        rxn <-rxn_df %>% pull(ID_old)
        if(i==1){
          rxn_name <- paste0("min_",rxn)
          p <- min_dataset_p %>% filter(metabolite == rxn_name) %>% pull(adjp)
          plot_df <- min_dataset[rxn_name]
        }else{
          rxn_name <- paste0("max_",rxn)
          p <- max_dataset_p %>% filter(metabolite == rxn_name) %>% pull(adjp)
          plot_df <- max_dataset[rxn_name]
        }
        plot_df$Cohort <- rep(c("before","pos"), each=49)
        plot_df$connect <- c(1:49,1:49)
        plots[[3*(i-1) + scale]] = ggplot(plot_df, aes_(x=as.name("Cohort"), y=as.name(rxn_name))) + geom_boxplot(aes(fill = Cohort))+geom_line(aes(group=connect), position = position_dodge())+geom_point(aes(fill=Cohort,group=connect),size=3,shape=21)+scale_fill_manual(labels = c(expression(italic("A.fumigatus -")), expression(italic("A.fumigatus +"))),values=c("#56B4E9","#E69F00"))+
  theme_pubr()  +annotate(geom = "text",label = paste("p =",formatC(p,3)),size=8,x=1.5,y=max(na.omit(plot_df[rxn_name]))*1.2)+
  theme(plot.title = element_text(size = 26, hjust = 0.5),  text = element_text(size = 24), legend.title = element_blank() , axis.text.x=element_blank(),axis.line = element_line(colour = 'black', size = 1),axis.ticks = element_line(colour = "black", size = 1),axis.ticks.length=unit(.2, "cm"),axis.text=element_text(size=20),
        legend.key.size = unit(2.3, "cm"),legend.text=element_text(size=24),legend.position = "none") +
  xlab(element_blank())+ylab(selected_rxn__final_name[scale])  }
}
plots <- list(plots[[1]]+ylim(-0.05, 0.05),plots[[4]]+ylab(""),plots[[2]]+ylim(-300, 0),plots[[5]]+ylab(""),plots[[3]],plots[[6]]+ylab(""))
g <- ggarrange(plotlist=plots,ncol=2,nrow = 3,align="hv")

ggsave("fig4D.pdf", width =8, height = 10,g)
```