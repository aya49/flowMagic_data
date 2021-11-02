# ----------------- generate csv file of groups ------------

generate_csv_groups_file<-function(list_groups,pop){
  vec<-unlist(list_groups)
  df<-data.frame(Samples=character(),Groups=character(),stringsAsFactors=FALSE)
  for (string in vec){
    str<-strsplit(string,":")[[1]]
    group<-str[1]
    sample<-str[2]
    new_row<-data.frame(sample,group)
    names(new_row)<-c("Samples","Groups")
    df<-rbind(df,new_row)
  }
  # export df to csv
  write.csv(df,file=sprintf("/home/rstudio/results/fcs_files_flowGroup/%s/df_groups.csv",pop),row.names = F)
  return(df)
}
