#-------------- import libraries ----------------
library(RchyOptimyx)
library(stringr)
library(flowTypeFilterC)
#----------- import csv file---

rchyoptimyx_df<-read.csv("/home/rstudio/data/rchyoptmix_data_test.csv",header = T)
flowTypeFilterC::encodePhenotype(as.character(str_remove(rchyoptimyx_df$Pheno_types[20],"_")), marker_names)
rchyoptimyx_results<-function(rchyoptimyx_df){
  pheno.codes<-rchyoptimyx_df$Pheno_types
  phenotypeScores<-rchyoptimyx_df$p_value
  inds<-which(phenotypeScores<0.05)
  selected_pvals<-phenotypeScores[inds]
  selected_pheno<-as.character(pheno.codes[inds])
  longest_selected_pheno<-selected_pheno[length(selected_pheno)]
  pheno.codes<-str_remove(pheno.codes,"_")
  pheno.codes<-str_remove(pheno.codes,"_")
  marker_names<-c(MarkerNames,filter_markers_bcells,filter_markers_blast,filter_markers_trans)
  marker_names<-marker_names[7:length(marker_names)]
  remove_markers<-c("HLADR","Time")
  string<-paste0(remove_markers,collapse = "|")
  inds<-grep(string,marker_names)
  marker_names<-marker_names[-inds]
  res<-RchyOptimyx(pheno.codes, -log10(phenotypeScores), startPhenotype=longest_selected_pheno,factorial(2),trimPaths = F)
  plot(res, phenotypeScores=-log10(phenotypeScores), phenotypeCodes=selected_pheno, marker.names=marker_names, ylab='-log10(Pvalue)')
}
rchyoptimyx_results(rchyoptimyx_df)