######### test CCP data ##########
f_example_data<-read.csv(file = "/home/rstudio/results/CLR_results/data/FCT_G509K_4R31.csv",check.names = F)
f_example_clrs<-read.csv(file = "/home/rstudio/results/CLR_results/clrs/FCT_G509K_4R31.csv",check.names = F)
nrow(f_example_clrs)
nrow(f_example_data)
color_col<-rep("a",nrow(f_example_data))
inds<-which(f_example_clrs[,5]==T)
color_col[inds]<-"popplus"
color_col[-inds]<-"popminus"
f_example_data<-cbind(f_example_data,color_col)
f_example_data$color_col<-as.factor(f_example_data$color_col)
ggplot(f_example_data,aes(x=f_example_data[,1],y=f_example_data[,2],colour=color_col)) + geom_point()



######### test CCP data ##########
f_example_data<-read.csv(file = "/home/rstudio/data/test_ccp_data/data/FCT_G001D_7F74.csv",check.names = F)
f_example_clrs<-read.csv(file = "/home/rstudio/data/test_ccp_data/clrs/FCT_G001D_7F74.csv",check.names = F)
nrow(f_example_clrs)
nrow(f_example_data)
color_col<-rep("a",nrow(f_example_data))
inds<-which(f_example_clrs$`gd+`==T)
color_col[-inds]<-"gdminus"
color_col[inds]<-"gdplus"
f_example_data<-cbind(f_example_data,color_col)
f_example_data$color_col<-as.factor(f_example_data$color_col)
ggplot(f_example_data,aes(x=`PE-CF594-A`,y=`PE-A`,colour=color_col)) + geom_point()
