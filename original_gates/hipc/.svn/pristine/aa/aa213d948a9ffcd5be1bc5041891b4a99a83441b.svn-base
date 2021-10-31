
#----- Gating margins----------------
margin_gating <- function(fs,gs){
  set.seed(40)
  margin.gates <- fsApply(fs,removeMargins,c("FSC-A","SSC-A"),return.gate=TRUE)
  rg <- lapply(1:nrow(margin.gates), function(x) return(rectangleGate(filterId="Margin", "FSC-A"=c(-1, margin.gates[x,1]),
                                                                      "SSC-A"=c(-1, margin.gates[x,2]))))
  names(rg) <-sampleNames(gs)
  node_margin<-add(gs, rg)#,"Margin")
  recompute(gs)
  fs.marg <- getData(gs,"Margin")
  return(list(gs=gs,fs.marg=fs.marg))
}




# #------- Compensation e transformation sul GatingSet  -------------------

comp_and_transform<- function(gs,comp){
  set.seed(40)
  gs<- compensate.flow(gs,comp = comp)
  recompute(gs)
  gs <-transform.flow(gs,remove.outliers=F,trans.chans = NULL)
  recompute(gs)
  fs.marg<-getData(gs,"Margin")
  return(list(gs=gs,fs.marg=fs.marg))
}


#------- Gating non marginal to cleaned ---------------------------
cleaning<-function(fs.marg,gs,n_cores,path.output){
  set.seed(40)
  print("Cleaning")
  clean.inds <- mclapply(1:length(fs.marg),function(x){
    message(paste("","Cleaning over time",sep="\n"))
    f<-fs.marg[[x]] 
    print(x) 
    segment <- 3000
    if (nrow(f)<5000)
      segment<-ceiling(nrow(f)/12) 
    
    cleaned <- cleaned.fcs(f = f, cleaning=T, transformed =T,id=identifier(f), segment=segment,comp.matrix=NULL, spill=F,PrintToConsole = F,make.NA=F,return.inds = F,return.all=T,
                           result.folder=  paste0(path.output,"/Cleaning/"))
    
    return(cleaned)
  },mc.cores = n_cores,mc.set.seed = F)
  names(clean.inds)<- sampleNames(gs)
  
  clean.clust <- mclapply(1:length(clean.inds), function(x){
    vec<-rep(0,nrow(fs.marg[[x]]))
    if(length(clean.inds[[x]]$ind)>0){
      vec[clean.inds[[x]]$ind]<-1
    }
    
    vec <- as.factor(vec)
    levels(vec) <- c("0", "1")
    return(vec)
  },mc.cores= n_cores)
  names(clean.clust)<-sampleNames(gs)
  add(gs,clean.clust,parent = "Margin", name = "Clean")
  recompute(gs)
  fs.clean <- getData(gs,"Clean_0")
  return(list(clean.inds=clean.inds,fs.clean=fs.clean,gs=gs))
}






# ---------------------------- Gating from cleaned cells to singlets --------------------------
gating_time_to_singlets<- function(fs.clean,gs,pre_processing,n_cores,scat.chans){
  print("gating_Time_to_singlets")
  set.seed(40)
  singlets_poly_list<-mclapply(1:length(gs),function(i){
    print(identifier(gs[[i]]@data[[i]]))
    theta0 <- pi/5.2 
    rot <- rotate.data(fs.clean[[i]], c(scat.chans['FSC-A'], scat.chans['FSC-H']), theta = theta0)$data
    
    rot.gate <- deGate(rot, channel = scat.chans['FSC-H'], use.upper = T, upper = F, alpha = 0.05)
    rot.gate_001 <- deGate(rot, channel = scat.chans['FSC-H'], use.upper = T, upper = F, alpha = 0.001)
    print(paste0("rot.gate:",rot.gate))
    print(paste0("rot.gate_001:",rot.gate_001))
    print(rot.gate-rot.gate_001)
    if((rot.gate-rot.gate_001)>60000){
      rot.gate<-rot.gate*1.20
      group_string<-paste0("group_1:",identifier(gs[[1]]@data[[i]]))
    }else if((rot.gate-rot.gate_001)>35000){
      rot.gate<-(rot.gate+rot.gate_001)/4
      group_string<-paste0("group_2:",identifier(gs[[1]]@data[[i]]))
    }else if((rot.gate-rot.gate_001)>10000){
      rot.gate<-(rot.gate+rot.gate_001)/2
      group_string<-paste0("group_3:",identifier(gs[[1]]@data[[i]]))
    }else{
      rot.gate<-rot.gate_001
      group_string<-paste0("group_4:",identifier(gs[[1]]@data[[i]]))
    }
    # show(autoplot(rot,"FSC-A","FSC-H"))
    # show(autoplot(rot,"FSC-H"))
    print(paste0("rot.gate:",rot.gate))
    if(rot.gate<(-70000)){
      rot.gate <- rot.gate/7
    }else if(rot.gate<(-10000)){
      rot.gate <- rot.gate/4
    }
    print(paste0("rot.gate:",rot.gate))
    singlets <- rotate.fd(flowDensity(rot, channels = c(scat.chans['FSC-A'], scat.chans['FSC-H']), position = c(NA, T), gates = c(NA, rot.gate)),angle = theta0)@filter
    singlets_pop <- rotate.fd(flowDensity(rot, channels = c(scat.chans['FSC-A'], scat.chans['FSC-H']), position = c(NA, T), gates = c(NA, rot.gate)),angle = theta0)
    
    sngl.poly <- polygonGate(filterId = "Singlets",.gate=singlets)
    return(list(sngl.poly=sngl.poly,singlets=singlets,groups=group_string,rot.gate=rot.gate,singlets_pop=singlets_pop))
  },mc.cores=n_cores,mc.set.seed = F)
  
  list_singlets_poly <- lapply(1:length(gs),function(i){
    singlets_poly_list[[i]]$sngl.poly
  })
  list_singlets <- lapply(1:length(gs),function(i){
    singlets_poly_list[[i]]$singlets
  })
  list_singlets_group <- lapply(1:length(gs),function(i){
    singlets_poly_list[[i]]$groups
  })
  list_rot.gate <- lapply(1:length(gs),function(i){
    singlets_poly_list[[i]]$rot.gate
  })
  list_singlets_pop <- lapply(1:length(gs),function(i){
    singlets_poly_list[[i]]$singlets_pop
  })
  names(list_singlets_poly)<-sampleNames(gs)
  if(pre_processing==T){
    nodeID0<-add(gs,list_singlets_poly,parent="Clean_0")
  }else{
    nodeID0<-add(gs,list_singlets_poly,parent="root")
  }
  recompute(gs)
  fs.sngl<-getData(gs,"Singlets")
  
  return(list(fs.sngl=fs.sngl,gs=gs,list_singlets=list_singlets,list_singlets_group=list_singlets_group,list_rot.gate=list_rot.gate,list_singlets_pop=list_singlets_pop))
}






# ------------------------ Gating from Singlets to Size and Beads----------------------


gating_singlets_to_size_beads<-function(fs.sngl,gs,n_cores,scat.chans){
  print("gating_Size_beads")
  # show(autoplot(fs.sngl[[1]],"SSC-A"))
  # show(autoplot(fs.sngl[[1]],"FSC-A","SSC-A"))
  set.seed(40)
  beads_size_poly_list<-mclapply(1:length(gs),function(i){
    print(identifier(gs[[i]]@data[[i]]))
    #theta0 <- pi/6.5 # defining rotation angle
    theta0 <- (1/12)*pi # defining rotation angle (15 degree rotation)
    rot <- rotate.data(fs.sngl[[i]], c(scat.chans['FSC-A'], scat.chans['SSC-A']), theta = theta0)$data
    #show(autoplot(rot,"FSC-A","SSC-A"))
    #show(autoplot(rot,"SSC-A"))
    all.cuts<-deGate(rot, channel = c(scat.chans['SSC-A']), all.cuts = T,tinypeak.removal = 0.01)
    all.peaks<-getPeaks(rot, channel = c(scat.chans['SSC-A']), tinypeak.removal = 0.01)
    print(paste0("all.cuts:",all.cuts))
    inds<-which(all.cuts<145000)
    all.cuts<-all.cuts[inds]
    max_cuts_gate<-max(all.cuts)
    if((max_cuts_gate<=50000)&&(max(all.peaks$Peaks)>190000)){ # beads peak is very high and cuts gate is too much low
      all.cuts<-deGate(rot, channel = c(scat.chans['SSC-A']), all.cuts = T,tinypeak.removal = 0.01)
      inds<-which(all.cuts<170000)
      all.cuts<-all.cuts[inds]
      max_cuts_gate<-max(all.cuts)
    }
    print(all.peaks)
    inds<-which(all.peaks$Peaks>110000)
    all.peaks$Peaks<-all.peaks$Peaks[inds]
    min_peak<-min(all.peaks$Peaks)
    print(paste0("min_peak:",min_peak))
    if(max_cuts_gate<90000){
      all.cuts<-deGate(rot, channel = c(scat.chans['SSC-A']), all.cuts = T,tinypeak.removal = 0.00000000001)
      #print(all.cuts)
      inds<-which(all.cuts<145000)
      all.cuts<-all.cuts[inds]
      max_cuts_gate<-max(all.cuts)
    }
    print(paste0("max_cuts_gate:",max_cuts_gate))
    if(max_cuts_gate<40000){
      max_cuts_gate<-deGate(rot, channel = c(scat.chans['SSC-A']), upper = T,tinypeak.removal = 0.0001,alpha = 0.2)
    }
    
    if(min_peak<120000){# the beads peak is too much in low position
      if(min_peak>110000){
        max_cuts_gate<-max_cuts_gate*0.80
      }else{
        max_cuts_gate<-max_cuts_gate*0.90
      }
    }
    
    print(paste0("max_cuts_gate:",max_cuts_gate))
    temp<- flowDensity(rot, c(scat.chans['FSC-A'], scat.chans['SSC-A']), position = c(NA, T), gates = c(NA,max_cuts_gate))
    #show(autoplot(getflowFrame(temp),"FSC-A","SSC-A"))
    #show(autoplot(getflowFrame(temp),"SSC-A"))
    rot.SSCA.gate<-deGate(temp, channel = c(scat.chans['SSC-A']), use.upper=T,upper=F,tinypeak.removal = 0.01,alpha = 0.001)
    if(max_cuts_gate>=135000){ # for the sample G5 plate 2
      rot.SSCA.gate<-rot.SSCA.gate*1.05
    }
    print(paste0("rot.SSCA.gate:",rot.SSCA.gate))
    beads<- rotate.fd(flowDensity(rot, c(scat.chans['FSC-A'], scat.chans['SSC-A']), position = c(NA, T), gates = c(NA,rot.SSCA.gate)),angle = theta0)
    #show(autoplot(getflowFrame(beads),"FSC-A","SSC-A"))
    #show(autoplot(getflowFrame(beads),"FSC-A"))
    #beads <- flowDensity(rot, channels = c(scat.chans['FSC-A'], scat.chans['SSC-A']), position = c(NA, T), gates = c(NA, rot.SSCA.gate*1.05))
    FSCAbead.gate<-deGate(beads, channel = c(scat.chans['SSC-A']), upper=T,alpha = 0.8)
    print(paste0("FSCAbead.gate:",FSCAbead.gate))
    if((FSCAbead.gate<70000)||(FSCAbead.gate>100000)){
      all.peaks<-getPeaks(beads, channel = c(scat.chans['FSC-A']),tinypeak.removal = 0.1)
      print(all.peaks)
      inds<-which(all.peaks$Peaks>40000)
      all.peaks$Peaks<-all.peaks$Peaks[inds]
      FSCAbead.gate<-min(all.peaks$Peaks)
      FSCAbead.gate<-FSCAbead.gate*1.15
      print(FSCAbead.gate)
    }
    
    print(paste0("FSCAbead.gate:",FSCAbead.gate))
    beads.poly<- polygonGate(filterId = "Beads",.gate=flowDensity(beads, channels = c(scat.chans['FSC-A'], scat.chans['SSC-A']), 
                                                                  position = c(F, NA), gates = c(FSCAbead.gate, NA), ellip.gate = T, scale = 0.96)@filter)
    beads<-flowDensity(beads, channels = c(scat.chans['FSC-A'], scat.chans['SSC-A']), position = c(F, NA), gates = c(FSCAbead.gate, NA), ellip.gate = T, scale = 0.96)
    # ------- finding the size population 
    
    #print(SSCA.gate)
    size<- rotate.fd(flowDensity(rot, c(scat.chans['FSC-A'], scat.chans['SSC-A']), position = c(NA,F), gates = c(NA, rot.SSCA.gate*0.98)),angle = theta0)
    #show(autoplot(getflowFrame(size),"FSC-A","SSC-A"))
    #size<- flowDensity(fs.sngl[[i]], channels = c(scat.chans['FSC-A'], scat.chans['SSC-A']), position = c(NA, F), gates = c(NA, SSCA.gate))
    theta0 <- pi/5.2 # definiamo angolo di rotazione
    rot <- rotate.data(getflowFrame(size), c(scat.chans['FSC-A'], scat.chans['SSC-A']), theta = theta0)$data
    #show(autoplot(rot,"FSC-A","SSC-A"))
    #show(autoplot(rot,"FSC-A"))
    rot.gate_all.cuts <- deGate(rot, channel = scat.chans['FSC-A'],upper = T, alpha = 0.05, all.cuts = T,tinypeak.removal=0.01)
    print(paste0("rot.gate_all.cuts:",rot.gate_all.cuts))
    rot.gate<-min(rot.gate_all.cuts)
    print(paste0("rot.gate_cuts:",rot.gate))
    all.peaks<-getPeaks(rot,channel = scat.chans['FSC-A'],tinypeak.removal = 0.01)
    print(all.peaks)
    inds<-which(all.peaks$Peaks>0)
    all.peaks$Peaks<-all.peaks$Peaks[inds]
    if((rot.gate<50000)&&(is.numeric(all.peaks$Peaks[2])==T)){
      rot.gate_1<-(all.peaks$Peaks[1]+all.peaks$Peaks[2])/2
      rot.gate<-mean(c(rot.gate_1,rot.gate))*0.95
    }
    left_peak<-min(all.peaks$Peaks)
    print(paste0("left_peak:",left_peak))
    if(left_peak>115000){ # left peak on the extreme right (maybe too many low cells in the file)
      group_string<-paste0("group_1:",identifier(gs[[i]]@data[[i]]))
      rot.gate <- deGate(temp, channel = scat.chans['FSC-A'],use.upper=F, upper = F)
      if(rot.gate>100000){# gate still incorrect
        rot.gate<-left_peak/2*0.90
      }
      print(paste0("temp_left_peak_100000:",rot.gate))
    }else if((left_peak>=47000)&&(left_peak<=115000)){ # it means that deGate and so getPeaks don't see the first peak or it is absent
      group_string<-paste0("group_2:",identifier(gs[[i]]@data[[i]]))
      rot.gate<-deGate(rot, channel = scat.chans['FSC-A'],use.upper=F,upper = F,tinypeak.removal = 0.01,alpha=0.9)
      print(paste0("rot.gate_1:",rot.gate))
      if(rot.gate>70000){
        temp<-flowDensity(rot, c(scat.chans['FSC-A'], scat.chans['SSC-A']), position = c(F, NA), gates = c(left_peak, NA))
        #show(autoplot(getflowFrame(temp),"FSC-A"))
        rot.gate<-deGate(temp, channel = scat.chans['FSC-A'],use.percentile = T,percentile = 0.5)
        rot.gate_1<-deGate(temp, channel = scat.chans['FSC-A'],upper=F)
        rot.gate<-mean(c(rot.gate,rot.gate_1))
        print(left_peak-rot.gate)
        print(paste0("temp:",rot.gate))
        if((left_peak-rot.gate)<20000){ # the peak and the gate are too close
          rot.gate<-deGate(temp, channel = scat.chans['FSC-A'],upper = F)*0.95
        }
      }
      if(rot.gate<31000){
        rot.gate<-(left_peak+rot.gate)/2
        if(left_peak>90000){
          rot.gate<-left_peak/2
        }
        
      }
    }else if((left_peak<47000)&&(length(all.peaks$Peaks)==1)){ # it means that getPeaks see only the first peak on the extreme left
      group_string<-paste0("group_3:",identifier(gs[[i]]@data[[i]]))
      rot.gate<-deGate(rot, channel = scat.chans['FSC-A'],after.peak = T,alpha = 0.9)
      if(rot.gate>65000){ # gate too much on the right
        rot.gate<-(rot.gate+left_peak)/2
      }
    }else{
      group_string<-paste0("group_4:",identifier(gs[[i]]@data[[i]]))
    }
    
    print(paste0("rot.gate:",rot.gate))
    if((rot.gate>85000)&&(left_peak<18000)){
      rot.gate<-left_peak*4
    }else if((rot.gate>85000)&&(left_peak>25000)){
      rot.gate<-left_peak*1.5
    }else if((rot.gate>85000)&&(left_peak<19000)){
      rot.gate<-left_peak*2
    }
    
    print(paste0("rot.gate:",rot.gate))
    size.temp<- rotate.fd(flowDensity(rot, c(scat.chans['FSC-A'], scat.chans['SSC-A']), position = c(T, NA), gates = c(rot.gate, NA)),angle = theta0)
    #print(size.temp)
    size.poly <- polygonGate(filterId = "All cells",.gate=size.temp@filter)
    return(list(beads.poly=beads.poly,size.poly=size.poly,beads=beads,groups=group_string,size.temp=size.temp))
  },mc.cores = n_cores,mc.set.seed = F)
  beads_poly_list<-lapply(1:length(gs),function(i){
    beads_size_poly_list[[i]]$beads.poly
  })
  size_poly_list<-lapply(1:length(gs),function(i){
    beads_size_poly_list[[i]]$size.poly
  })
  
  list_beads<-lapply(1:length(gs),function(i){
    beads_size_poly_list[[i]]$beads
  })
  list_groups_size<-lapply(1:length(gs),function(i){
    beads_size_poly_list[[i]]$groups
  })
  list_size<-lapply(1:length(gs),function(i){
    beads_size_poly_list[[i]]$size.temp
  })
  names(beads_poly_list)<-sampleNames(gs)
  nodeID0<-add(gs,beads_poly_list,parent="Singlets")
  recompute(gs)
  names(size_poly_list)<-sampleNames(gs)
  nodeID0<-add(gs,size_poly_list,parent="Singlets")
  recompute(gs)
  fs.size<-getData(gs,"All cells")
  return(list(fs.size=fs.size,gs=gs,beads_poly_list=beads_poly_list,size_poly_list=size_poly_list,list_beads=list_beads,list_groups_size=list_groups_size,list_size=list_size))
}
# ------------------------ Gating from size to Live cells ----------------------

gating_size_to_live<-function(fs.size,gs,n_cores,fluorochrome.chans,scat.chans){
  print("gating_live")
  set.seed(40)
  live_poly_list<-mclapply(1:length(gs),function(i){
    print(identifier(gs[[i]]@data[[i]]))
    # in order to obtain an oblique gates,we need to rotate the data,
    # but we cannot rotate using big angles,because the data are scaled 
    # in a different way along the x axis and y axis,
    # I use the min.max option of the rotate.data function 
    #to find the right angle rotation, and 
    # I play around this angle(only a little bit,the numbers must be very small)
    #show(autoplot(fs.size[[1]],"APC-eF780-A","SSC-A"))
    output <- rotate.data(fs.size[[i]], c(fluorochrome.chans["Live"], scat.chans['SSC-A']), min.max = T)
    theta0 <- output$theta
    theta0<-theta0/8
    rot <- rotate.data(fs.size[[i]], c(fluorochrome.chans["Live"], scat.chans['SSC-A']), theta = theta0)$data
    # show(autoplot(rot,"APC-eF780-A","SSC-A"))
    # show(autoplot(rot,"APC-eF780-A"))
    rot.live.gate_all.cuts <- deGate(rot,"APC-eF780-A",all.cuts=T)
    all.peaks <- getPeaks(rot,"APC-eF780-A",tinypeak.removal = 0.1)
    rot.live.gate_all.peaks<-all.peaks$Peaks
    print(paste0("all.peaks.rot:",rot.live.gate_all.peaks))
    print(paste0("all.cuts.rot:",rot.live.gate_all.cuts))
    max_peak_rot<-max(rot.live.gate_all.peaks)
    print(paste0("max_peak_rot:",max_peak_rot))
    if(max_peak_rot>=2.81){ # the max peak is on the right,we put the gate on the left of this peak
      group_string<-paste0("group_1:",identifier(gs[[1]]@data[[i]]))
      rot.temp.live.gate<-max_peak_rot
      temp<-flowDensity(rot, c(fluorochrome.chans['Live'], scat.chans['SSC-A']), position = c(F, NA), gates = c(max_peak_rot, NA))
      rot.live.gate<-deGate(temp, channel = fluorochrome.chans['Live'], use.upper=T, upper = T, alpha = 0.04)
      rot.live.gate_2<-max_peak_rot*0.87
      if(rot.live.gate_2<3){
        rot.live.gate<-rot.live.gate_2
      }else{
        rot.live.gate<-rot.live.gate
      }
      # show(autoplot(getflowFrame(temp),"APC-eF780-A","SSC-A"))
      # show(autoplot(getflowFrame(temp),"APC-eF780-A"))
      print(rot.live.gate)
      if(rot.live.gate>=2.9){
        rot.live.gate<-deGate(temp, channel = fluorochrome.chans['Live'], use.upper=T,upper = T, alpha = 0.01,tinypeak.removal = 0.5)
      }
      print(rot.live.gate)
      if((rot.live.gate>=3.0)||(rot.live.gate<2)){
        rot.live.gate<-max(rot.live.gate_all.cuts)
        inds<-which(all.peaks$Peaks<rot.live.gate)
        if(length(inds)!=0){
          all.peaks$Peaks<-all.peaks$Peaks[inds]
          peaks_near_the_cut<-max(all.peaks$Peaks)
          diff<-rot.live.gate-peaks_near_the_cut
          if(diff<0.2){
            rot.live.gate<-rot.live.gate*1.04
          }
        }  
      }
      print(rot.live.gate)
      if(rot.live.gate>=2.9){
        all.peaks <- getPeaks(temp,"APC-eF780-A",tinypeak.removal = 0.001)
        print(all.peaks)
        inds<-which(all.peaks$Peaks<3)
        all.peaks.h<-all.peaks$P.h[inds]
        max.h.peak<-max(all.peaks.h)
        if(max.h.peak<0.2){
          rot.live.gate<-deGate(temp, channel = fluorochrome.chans['Live'], upper = F, alpha = 0.3,tinypeak.removal = 0.3)
        }
      }else if((rot.live.gate<=1.85)&&(length(rot.live.gate_all.peaks)==1)){
        rot.live.gate<-max_peak_rot*0.90
      }
      print(rot.live.gate)
      
    }else{ # the max peak is more on the left,we put the gate on the right of this peak
      group_string<-paste0("group_2:",identifier(gs[[1]]@data[[i]]))
      rot.live.gate<-deGate(rot, channel = fluorochrome.chans['Live'], use.upper=T, upper = T, alpha = 0.01)
      print(paste0("rot.live.gate_0.01:",rot.live.gate))
      rot.live.gate_2<-max_peak_rot*1.10
      if((rot.live.gate_2-max_peak_rot)<0.3){ # too much near
        rot.live.gate<-max(rot.live.gate_2,rot.live.gate)
      }else{ # enough far
        rot.live.gate<-min(rot.live.gate_2,rot.live.gate)
      }
      #----------------------------
      # gate too much on the right
      if(rot.live.gate>3.2){
        rot.live.gate<-deGate(rot, channel = fluorochrome.chans['Live'], use.upper=T, upper = T,alpha = 0.1)
        print(paste0("rot.live.gate_0.1:",rot.live.gate))
        if((max_peak_rot>2.3)&&(max_peak_rot<2.7)){
          rot.live.gate<-deGate(rot, channel = fluorochrome.chans['Live'], use.upper=T, upper = T, alpha = 0.9)
          print(paste0("rot.live.gate_0.9:",rot.live.gate))
        }
        if(max_peak_rot<2.3){
          rot.live.gate<-deGate(rot, channel = fluorochrome.chans['Live'], use.upper=T, upper = T, alpha = 0.9)
          rot.live.gate<-rot.live.gate*1.10
          print(paste0("rot.live.gate_0.9_+:",rot.live.gate))
        }
      }
      #--------------------------- 
      # gate still too much on the right
      if((rot.live.gate>3.2)&&(max_peak_rot<2.7)){
        inds<-which(rot.live.gate_all.cuts<3.1)
        rot.live.gate_all.cuts<-rot.live.gate_all.cuts[inds]
        rot.live.gate<-max(rot.live.gate_all.cuts)
        print(paste0("rot.live.gate_cuts:",rot.live.gate))
      }
      #-------------------------------- 
      # gate too much on the left
      if((rot.live.gate<1.65)&&(max_peak_rot>=1)){
        rot.live.gate<-deGate(rot, channel = fluorochrome.chans['Live'], use.upper=T, upper = T, alpha = 0.2,tinypeak.removal = 0.3)
        print(paste0("rot.live.gate_0.2:",rot.live.gate))
      } else if((max_peak_rot<1)&&(rot.live.gate<1.65)){
        rot.live.gate.0.18<-deGate(rot, channel = fluorochrome.chans['Live'], use.upper=T, upper = T,alpha = 0.18)
        rot.live.gate.0.01<-deGate(rot, channel = fluorochrome.chans['Live'], use.upper=T, upper = T, alpha = 0.01)
        rot.live.gate <- mean(c(rot.live.gate.0.18,rot.live.gate.0.01))
        print(paste0("rot.live.gate_mean:",rot.live.gate))
      }else if((rot.live.gate<2.05)&&(max_peak_rot<=2.1)){
        rot.live.gate<-rot.live.gate*1.10
      }
    }
    
    print(paste0("rot.live.gate:",rot.live.gate))
    live<-rotate.fd(flowDensity(rot, c(fluorochrome.chans['Live'], scat.chans['SSC-A']), position = c(F, NA), gates = c(rot.live.gate, NA)),angle = theta0)
    
    
    #show(autoplot(getflowFrame(not_live),"APC-eF780-A","SSC-A"))
    #show(autoplot(getflowFrame(live),"APC-eF780-A","SSC-A"))
    
    #live.filter <-flowDensity(fs.size[[i]],channels =c(fluorochrome.chans['Live'], scat.chans['SSC-A']), position=c(F,NA),gates=c(live.gate,NA))@filter
    live.poly <- polygonGate(filterId = "Live cells",.gate  =live@filter)
    return(list(live.poly=live.poly,groups=group_string,live=live))
  },mc.cores = n_cores,mc.set.seed = F)
  lista_poly_live<-lapply(1:length(gs),function(i){
    live_poly_list[[i]]$live.poly
  })
  lista_groups_live<-lapply(1:length(gs),function(i){
    live_poly_list[[i]]$groups
  })
  list_live<-lapply(1:length(gs),function(i){
    live_poly_list[[i]]$live
  })
  names(lista_poly_live)<-sampleNames(gs)
  nodeID0<-add(gs,lista_poly_live,parent="All cells")
  recompute(gs)
  fs.live <- getData(gs, "Live cells")
  return(list(fs.live=fs.live,gs=gs,lista_poly_live=lista_poly_live,lista_groups_live=lista_groups_live,list_live=list_live))
}
# ------------------------ Gating from Live to Non Granulocytes----------------------
gating_live_to_non_Granulocytes<- function(fs.live,gs,n_cores,fluorochrome.chans){
  print("gating_to_non_granulocytes")
  set.seed(40)
  list_poly_non_gran<-mclapply(1:length(gs),function(i){
    # show(autoplot(fs.live[[i]],"BV711-A","V500-A"))
    # show(autoplot(fs.live[[i]],"BV711-A"))
    cd66.gate <- deGate(fs.live[[i]],fluorochrome.chans['CD66'],upper=T,tinypeak.removal = 0.1)
    print(cd66.gate)
    if(cd66.gate<1.5){ # gate too on left
      cd66.gate<-cd66.gate*1.50
    }
    if(cd66.gate>3){ # too on the right
      cd66.gate<-deGate(fs.live[[i]],fluorochrome.chans['CD66'],upper=T) # tinypeak removal removed. 
    }
    print(paste0("CD66:",cd66.gate))
    group_string<-paste0("group_1:",identifier(gs[[1]]@data[[i]]))
    non_gran<- flowDensity(fs.live[[i]], channels = c(fluorochrome.chans['CD66'], fluorochrome.chans['CD14']), position = c(F, NA), gates = c(cd66.gate, NA))
    gran<- flowDensity(fs.live[[i]], channels = c(fluorochrome.chans['CD66'], fluorochrome.chans['CD14']), position = c(T, NA), gates = c(cd66.gate, NA))
    non.gran.poly <- polygonGate(filterId = "Non granulocytes",.gate=non_gran@filter)
    gran.poly <- polygonGate(filterId = "Granulocytes",.gate  =gran@filter)
    return(list(non.gran.poly=non.gran.poly,gran.poly=gran.poly,groups=group_string,gran=gran,non_gran=non_gran))
  },mc.cores=n_cores)
  #--------------------------------------------
  lista_poly_non_gran<-lapply(1:length(gs),function(i){
    list_poly_non_gran[[i]]$non.gran.poly
  })
  lista_poly_gran<-lapply(1:length(gs),function(i){
    list_poly_non_gran[[i]]$gran.poly
  })
  lista_groups_gran<-lapply(1:length(gs),function(i){
    list_poly_non_gran[[i]]$groups
  })
  list_gran<-lapply(1:length(gs),function(i){
    list_poly_non_gran[[i]]$gran
  })
  list_non_gran<-lapply(1:length(gs),function(i){
    list_poly_non_gran[[i]]$non_gran
  })
  names(lista_poly_non_gran)<-sampleNames(gs)
  nodeID0<-add(gs,lista_poly_non_gran,parent="Live cells")
  names(lista_poly_gran)<-sampleNames(gs)
  nodeID0<-add(gs,lista_poly_gran,parent="Live cells")
  recompute(gs) # modifichiamo effettivamente il gatingSet in base al nuovo albero
  fs.gran <- getData(gs, "Granulocytes")
  fs.non_gran <-getData(gs,"Non granulocytes")
  return(list(fs.gran=fs.gran,gs=gs,fs.non_gran=fs.non_gran,lista_poly_gran=lista_poly_gran,lista_poly_non_gran=lista_poly_non_gran,lista_groups_gran=lista_groups_gran,list_non_gran=list_non_gran,list_gran=list_gran))
  
  #  show(autoplot(fs.live[[1]],"BV711-A")), show(autoplot(fs.live[[1]],"V500-A"))
  #  show(autoplot(fs.live[[1]],"BV711-A","V500-A"))
}




# ------------------------ Gating from non-Granulocytes to Bcells(CD19+)----------------------
# show(autoplot(fs.non_gran[[2]],"APC-A","SSC-A"))
# show(autoplot(fs.non_gran[[33]],"APC-A","BV605-A"))
# show(autoplot(fs.non_gran[[1]],"APC-A"))
# show(autoplot(fs.non_gran[[3]],"SSC-A"))
gating_nonGran_to_Bcells<- function(fs.non_gran,gs,fluorochrome.chans,scat.chans){
  print("gating_Nongran_to_Bcells")
  set.seed(40)
  list_Bcells_poly<-lapply(1:length(gs), function(i){
    #--- rotate data
    output <- rotate.data(fs.non_gran[[i]], c(fluorochrome.chans["CD19"],scat.chans["SSC-A"]), min.max = T)
    theta0 <- output$theta
    theta0<-theta0/2
    rot <- rotate.data(fs.non_gran[[i]],  c(fluorochrome.chans["CD19"],scat.chans["SSC-A"]), theta = theta0)$data
    #show(autoplot(rot,"APC-A","SSC-A"))
    # show(autoplot(rot,"APC-A"))
    #--- find the gate in the rotated data
    cd19.gate <- deGate(rot, fluorochrome.chans["CD19"],upper=T,tinypeak.removal=0.01,alpha = 0.1)
    print(paste0("cd19.gate:",cd19.gate))
    all.peaks<-getPeaks(rot,fluorochrome.chans["CD19"],tinypeak.removal=0.01)
    print(all.peaks)
    if(max(all.peaks$Peaks)<1.8){
      cd19.gate <- deGate(rot, fluorochrome.chans["CD19"], use.upper=T, upper=T, tinypeak.removal=0.01,alpha = 0.1) 
    }
    if(cd19.gate<1.6 || cd19.gate>2.5){
      if(max(all.peaks$Peaks)<2){
        cd19.gate <- deGate(rot, fluorochrome.chans["CD19"], use.upper=T, upper=T, tinypeak.removal=0.01,alpha = 0.8) 
      }else if(max(all.peaks$Peaks)>2){
        cd19.gate<-max(all.peaks$Peaks)*0.85
      }
      print(paste0("cd19.gate:",cd19.gate))
    }
    if(cd19.gate>2.5 && max(all.peaks$Peaks)<1.8){
      cd19.gate<-(max(all.peaks$Peaks)+cd19.gate)/2
    }
    print(paste0("cd19.gate:",cd19.gate))
    Bcells<-rotate.fd(flowDensity(rot,c(fluorochrome.chans["CD19"],scat.chans["SSC-A"]),position = c(T,NA),gates=c(cd19.gate,NA)),angle=theta0)
    NotBcells<-rotate.fd(flowDensity(rot,c(fluorochrome.chans["CD19"],scat.chans["SSC-A"]),position = c(F,NA),gates=c(cd19.gate,NA)),angle=theta0)
    
    Bcells_poly<- polygonGate(filterId = "CD19+ B cells",.gate=Bcells@filter)
    return(list(Bcells_poly=Bcells_poly,Bcells=Bcells,NotBcells=NotBcells))
  })
  #----------------------------------------------------------
  list_poly_Bcells<-lapply(1:length(gs),function(i){
    list_Bcells_poly[[i]]$Bcells_poly
  })
  list_groups_Bcells<-lapply(1:length(gs),function(i){
    list_Bcells_poly[[i]]$groups
  })
  list_Bcells<-lapply(1:length(gs),function(i){
    list_Bcells_poly[[i]]$Bcells
  })
  list_NotBcells<-lapply(1:length(gs),function(i){
    list_Bcells_poly[[i]]$NotBcells
  })
  names(list_poly_Bcells)<-sampleNames(gs)
  nodeID0<-add(gs,list_poly_Bcells,parent="Non granulocytes")
  recompute(gs) 
  fs.bcells<-getData(gs, "CD19+ B cells")
  
  return(list(gs=gs,fs.bcells=fs.bcells,list_poly_Bcells=list_poly_Bcells,list_groups_Bcells=list_groups_Bcells,
              list_Bcells=list_Bcells,list_NotBcells=list_NotBcells))
}


# ------------------------ Gating from Bcells to PC----------------------
# show(autoplot(fs.19[[3]],"PE-Cy5-A","V450-A"))
# show(autoplot(fs.19[[1]],"PE-Cy5-A"))
# show(autoplot(fs.19[[2]],"V450-A"))
gating_bcells_to_pc<-function(fs.bcells,gs,fluorochrome.chans){
  print("gating_pc")
  set.seed(40)
  list_pc_poly<-lapply(1:length(gs),function(i){
    print(identifier(gs[[i]]@data[[i]]))
    cd138.gate_0.1 <- deGate(fs.bcells[[i]],fluorochrome.chans["CD138"],upper=T,use.upper=T,alpha = 0.1) # gate on th extreme high part of the graph
    cd138.gate_0.9 <- deGate(fs.bcells[[i]],fluorochrome.chans["CD138"],upper=T,use.upper=T,alpha = 0.9,tinypeak.removal = 0.1) # gate more down (but above the bigger population)
    if(cd138.gate_0.9>2.8){
      cd138.gate_0.9 <- deGate(fs.bcells[[i]],fluorochrome.chans["CD138"],upper=T,alpha = 0.9,tinypeak.removal = 0.1) # gate more down (but above the bigger population)
    }
    print(paste0("cd138.gate_0.1:",cd138.gate_0.1))
    print(paste0("cd138.gate_0.9:",cd138.gate_0.9))
    cd138.gate<-(cd138.gate_0.1+cd138.gate_0.9)/2
    print(paste0("CD138.gate:",cd138.gate))
    if((cd138.gate>2.4)&&(cd138.gate_0.9<2)){ # gate too much high
      group_string<-paste0("group_1:",identifier(gs[[1]]@data[[i]]))
      
      cd138.gate<-cd138.gate_0.9
      # if still too much high
      if(cd138.gate>2.4){
        cd138.gate <- deGate(fs.bcells[[i]],fluorochrome.chans["CD138"],upper=T,alpha = 0.9)
        if(cd138.gate<1){
          cd138.gate<-cd138.gate_0.9
        }
      }else if((cd138.gate<2.2)&&(cd138.gate_0.9>1.8)){ # or too much low
        cd138.gate<-(cd138.gate_0.1+cd138.gate_0.9)/2
      }
    }else if((cd138.gate<1.1)&&(cd138.gate_0.9>=0)){ # gate too much low
      cd138.gate<-cd138.gate_0.1
      group_string<-paste0("group_2:",identifier(gs[[1]]@data[[i]]))
      
    }else if((cd138.gate_0.1<2.1)&&(cd138.gate_0.9>=0)){ #  cd138.gate_0.1 is ok to be used as gate
      cd138.gate<-cd138.gate_0.1
      group_string<-paste0("group_3:",identifier(gs[[1]]@data[[i]]))
      print(paste0("mean_CD138_values:",mean(fs.bcells[[i]]@exprs[,fluorochrome.chans["CD138"]])))
      if((mean(fs.bcells[[i]]@exprs[,fluorochrome.chans["CD138"]])<0.5)&&(cd138.gate_0.1>1.3)){
        cd138.gate<-cd138.gate_0.9
      }
    }else if((cd138.gate>2.7)&&(cd138.gate_0.9>2.6)){  #  cd138.gate_0.9 is ok to be used as gate
      cd138.gate<-cd138.gate_0.9
      group_string<-paste0("group_4:",identifier(gs[[1]]@data[[i]]))
    }else{
      group_string<-paste0("group_5:",identifier(gs[[1]]@data[[i]]))
    }
    
    print(paste0("CD138.gate:",cd138.gate))
    cd38.gate  <- deGate(fs.bcells[[i]], fluorochrome.chans["CD38"],upper=F,alpha = 0.01)
    if(cd38.gate>2){
      cd38.gate  <- deGate(fs.bcells[[i]], fluorochrome.chans["CD38"],upper=F,use.upper=T,alpha = 0.8)
    }
    print(paste0("CD38:",cd38.gate))
    pc <- flowDensity(fs.bcells[[i]], channels= c(fluorochrome.chans["CD38"],fluorochrome.chans["CD138"]),position = c(T,T),gates=c(cd38.gate,cd138.gate))
    PC_poly<-polygonGate(filterId = "plasma cells",.gate=pc@filter)
    return(list(PC_poly=PC_poly,groups=group_string,pc=pc))
  })
  #--------------------------
  list_poly_pc<-lapply(1:length(gs),function(i){
    list_pc_poly[[i]]$PC_poly
  })
  list_groups_pc<-lapply(1:length(gs),function(i){
    list_pc_poly[[i]]$groups
  })
  list_pc<-lapply(1:length(gs),function(i){
    list_pc_poly[[i]]$pc
  })
  names(list_poly_pc)<-sampleNames(gs)
  nodeID0<-add(gs,list_poly_pc,parent="CD19+ B cells")
  recompute(gs)
  fs.pc<-getData(gs, "plasma cells")
  return(list(gs=gs,fs.pc=fs.pc,list_poly_pc=list_poly_pc,list_groups_pc=list_groups_pc,list_pc=list_pc))
}



# ------------------------Gating from Bcells to PB----------------------
#show(autoplot(fs.19[[1]],"PE-Cy5-A","PE-Cy7-A"))
#show(autoplot(fs.19[[3]],"PE-Cy5-A"))
#show(autoplot(fs.19[[1]],"PE-Cy7-A"))
gating_bcells_to_PB<-function(fs.bcells,gs,fluorochrome.chans){
  print("gating_PB")
  set.seed(40)
  list_PB_poly<-lapply(1:length(gs), function(i){
    cd38.gate <-deGate(fs.bcells[[i]],fluorochrome.chans["CD38"],upper = F,alpha = 0.01)
    if(cd38.gate>2){
      cd38.gate  <- deGate(fs.bcells[[i]], fluorochrome.chans["CD38"],upper=F,use.upper=T,alpha = 0.8)
    }
    group_string<-paste0("group_1:",identifier(gs[[1]]@data[[i]]))
    print(paste0("cd38.gate:",cd38.gate))
    cd27.gate_0.1<- deGate(fs.bcells[[i]], fluorochrome.chans["CD27"],use.upper=T,upper=T, alpha = 0.1,tinypeak.removal = 0.5)*1.05
    cd27.gate_0.9 <- deGate(fs.bcells[[i]], fluorochrome.chans["CD27"],use.upper=T,upper=T, alpha = 0.9,tinypeak.removal = 0.5)*1.05
    print(paste0("cd27.gate_0.1:",cd27.gate_0.1))
    print(paste0("cd27.gate_0.9:",cd27.gate_0.9))
    cd27.gate<-(cd27.gate_0.1+cd27.gate_0.9)/2
    if(cd27.gate>2.3){
      cd27.gate<-cd27.gate*0.95
    }
    print(paste0("cd27.gate:",cd27.gate))
    PB <- flowDensity(fs.bcells[[i]],channels=c(fluorochrome.chans["CD38"],fluorochrome.chans["CD27"]),position = c(T,T),gates=c(cd38.gate,cd27.gate))
    PB.poly<-polygonGate(filterId = "Plasmablasts",.gate=PB@filter)
    return(list(PB.poly=PB.poly,groups=group_string,PB=PB))
  })
  #-----------------------------------------------
  list_poly_PB<-lapply(1:length(gs),function(i){
    list_PB_poly[[i]]$PB.poly
  })
  list_groups_PB<-lapply(1:length(gs),function(i){
    list_PB_poly[[i]]$groups
  })
  list_PB<-lapply(1:length(gs),function(i){
    list_PB_poly[[i]]$PB
  })
  names(list_poly_PB)<-sampleNames(gs)
  nodeID0<-add(gs,list_poly_PB,parent="CD19+ B cells")
  recompute(gs)
  fs.PB<-getData(gs,"Plasmablasts")
  return(list(gs=gs,fs.PB=fs.PB,list_poly_PB=list_poly_PB,list_groups_PB=list_groups_PB,list_PB=list_PB))
}



# ------------------------Gating from Bcells to CD19/20----------------------
# show(autoplot(fs.19[[6]],"APC-A")), show(autoplot(fs.19[[3]],"BV786-A"))
# show(autoplot(fs.19[[1]],"APC-A","BV786-A"))
#show(autoplot(fs.19[[2]],"BV786-A"))

gating_bcells_to_CD20p<-function(fs.bcells,gs,fluorochrome.chans){
  print("gating_CD19/20")
  set.seed(40)
  list_CD19_20_poly<-lapply(1:length(gs), function(i){
    print(identifier(gs[[i]]@data[[i]]))
    # show(autoplot(fs.bcells[[i]],"APC-A","BV786-A"))
    # show(autoplot(fs.bcells[[i]],"BV786-A"))
    cd19.gate <-deGate(fs.bcells[[i]],fluorochrome.chans["CD19"],upper = T,use.upper=T,alpha=0.00005) # al posto di 0.5, ho messo un intensità di picco di 0.00005 in modo che il treshold sia più a destra(considero picchi di intensità più deboli come veri e propri picchi)
    cd20.gate <- deGate(fs.bcells[[i]], fluorochrome.chans["CD20"],upper=F,alpha = 0.1)
    print(cd20.gate)
    if(cd20.gate>3){
      cd20.gate <- deGate(fs.bcells[[i]], fluorochrome.chans["CD20"],use.upper=F,upper=F,alpha = 0.1,tinypeak.removal = 0.5)
    }
    print(cd20.gate)
    cd20.pos<- flowDensity(fs.bcells[[i]], channels= c(fluorochrome.chans["CD19"],fluorochrome.chans["CD20"]),position = c(F,T),gates=c(cd19.gate,cd20.gate))
    cd20.neg<- flowDensity(fs.bcells[[i]], channels= c(fluorochrome.chans["CD19"],fluorochrome.chans["CD20"]),position = c(F,F),gates=c(cd19.gate,cd20.gate))
    cd20.pos.poly<-polygonGate(filterId = "CD19+CD20+",.gate=cd20.pos@filter)
    cd20.neg.poly<-polygonGate(filterId = "CD19+CD20-",.gate=cd20.neg@filter)
    group_string<-paste0("group_1:",identifier(gs[[1]]@data[[i]]))
    
    
    return(list(cd20.pos.poly=cd20.pos.poly, cd20.neg.poly=cd20.neg.poly,groups=group_string,cd20.pos=cd20.pos,cd20.neg=cd20.neg))
  })
  #------------------------------------
  
  list_CD20_pos_poly<-lapply(1:length(gs),function(i){
    list_CD19_20_poly[[i]]$cd20.pos.poly
  })
  list_CD20_neg_poly<-lapply(1:length(gs),function(i){
    list_CD19_20_poly[[i]]$cd20.neg.poly
  })
  list_groups_CD20<-lapply(1:length(gs),function(i){
    list_CD19_20_poly[[i]]$groups
  })
  list_cd20.pos<-lapply(1:length(gs),function(i){
    list_CD19_20_poly[[i]]$cd20.pos
  })
  list_cd20.neg<-lapply(1:length(gs),function(i){
    list_CD19_20_poly[[i]]$cd20.neg
  })
  names(list_CD20_pos_poly)<-sampleNames(gs)
  names(list_CD20_neg_poly)<-sampleNames(gs)
  nodeID0<-add(gs,list_CD20_pos_poly,parent="CD19+ B cells")
  nodeID0<-add(gs,list_CD20_neg_poly,parent="CD19+ B cells")
  recompute(gs)
  fs.CD20pos<-getData(gs,"CD19+CD20+")
  fs.CD20neg<-getData(gs,"CD19+CD20-")
  
  return(list(fs.CD20pos=fs.CD20pos, gs=gs,fs.CD20neg=fs.CD20neg,list_CD20_pos_poly=list_CD20_pos_poly,
              list_CD20_neg_poly=list_CD20_neg_poly,list_groups_CD20=list_groups_CD20,list_cd20.pos=list_cd20.pos,list_cd20.neg=list_cd20.neg))
}


# ------------------------Gating from Bcells to IgM/IgD----------------------

gating_bcells_to_igm_igd<-function(fs.bcells,gs,fluorochrome.chans){
  print("gating_igm_igd")
  set.seed(40)
  list_All_quads_poly<-lapply(1:length(gs), function(i){
    print(identifier(gs[[i]]@data[[i]]))
    # show(autoplot(fs.bcells[[i]],"PE-CF594-A","FITC-A"))
    # show(autoplot(fs.bcells[[i]],"PE-CF594-A"))
    igm.gate<-deGate(fs.bcells[[i]],fluorochrome.chans["IgM"],upper=F,tinypeak.removal=0.1)
    igd.gate<-deGate(fs.bcells[[i]],fluorochrome.chans["IgD"],upper=F,tinypeak.removal = 0.6,alpha=0.8)
    print(paste0("igd_gate:",igd.gate))
    print(paste0("igm.gate:",igm.gate))
    all.peaks<-getPeaks(fs.bcells[[i]],fluorochrome.chans["IgD"],tinypeak.removal = 0.6)
    print(all.peaks)
    if(igd.gate>2.2 && max(all.peaks$Peaks)>1.8){
      all.peaks<-getPeaks(fs.bcells[[i]],fluorochrome.chans["IgD"],tinypeak.removal = 0.6)
      inds<-which(all.peaks$Peaks>2)
      all.peaks$Peaks<-all.peaks$Peaks[inds]
      igd.gate<-min(all.peaks$Peaks)*0.80
    }
    if(max(all.peaks$Peaks)<1.8){
        igd.gate<-deGate(fs.bcells[[i]],fluorochrome.chans["IgD"],upper=F,use.upper=F,alpha=0.8,tinypeak.removal = 0.6)
        print(igd.gate)
        if(igd.gate>0 && min(all.peaks$Peaks)<0.20 && igd.gate<0.25){
          igd.gate<-(-igd.gate)
        }
    }
    if(igd.gate<1 && max(all.peaks$Peaks)>1.8){
        all.cuts<-deGate(fs.bcells[[i]], channel = c(fluorochrome.chans["IgD"]), all.cuts = T,tinypeak.removal = 0.01)
        print(all.cuts)
        igd.gate<-max(all.cuts)
    }
    print(paste0("igd.gate:",igd.gate))
    temp <- flowDensity(fs.bcells[[i]], channels=c(fluorochrome.chans["IgD"],fluorochrome.chans["IgM"]),position=c(T,NA),gates=c(igd.gate,NA))
    igm.gate <- deGate(temp@flow.frame, fluorochrome.chans["IgM"],use.upper=T,upper=F,alpha=0.5)
    if(igm.gate<1.5){
      all.peaks<-getPeaks(fs.bcells[[i]],fluorochrome.chans["IgM"])
      max_peak<-max(all.peaks$Peaks)
      igm.gate<-(max_peak+igm.gate)/2
    }
    print(paste0("igm.gate:",igm.gate))
    quad.1 <- flowDensity(fs.bcells[[i]], channels=c(fluorochrome.chans["IgD"],fluorochrome.chans["IgM"]),position=c(F,F),gates=c(igd.gate,igm.gate))
    quad.2 <- flowDensity(fs.bcells[[i]],channels=c(fluorochrome.chans["IgD"],fluorochrome.chans["IgM"]),position=c(T,F), gates=quad.1@gates)
    quad.3 <- flowDensity(fs.bcells[[i]], channels=c(fluorochrome.chans["IgD"],fluorochrome.chans["IgM"]),position=c(F,T), gates=quad.1@gates)
    quad.4 <- flowDensity(fs.bcells[[i]], channels=c(fluorochrome.chans["IgD"],fluorochrome.chans["IgM"]),position=c(T,T), gates=quad.1@gates)
    igm <- flowDensity(fs.bcells[[i]], channels=c(fluorochrome.chans["IgD"],fluorochrome.chans["IgM"]),position=c(NA,T), gates=quad.1@gates)
    quad.1.poly<-polygonGate(filterId = "IGD-IGM- B cells",.gate=quad.1@filter)
    quad.2.poly<-polygonGate(filterId = "IGD+IGM- B cells",.gate=quad.2@filter)
    quad.3.poly<-polygonGate(filterId = "IGD-IGM+ B cells",.gate=quad.3@filter)
    quad.4.poly<-polygonGate(filterId = "IGD+IGM+ B cells",.gate=quad.4@filter)
    igm.poly<-polygonGate(filterId = "IGM",.gate=igm@filter)
    return(list(quad.1.poly=quad.1.poly,quad.2.poly=quad.2.poly,quad.3.poly=quad.3.poly, quad.4.poly=quad.4.poly, 
                igm.poly=igm.poly,quad.1=quad.1,quad.1=quad.1,quad.2=quad.2,quad.3=quad.3,quad.4=quad.4))
    
  })
  #-----------------------------------------------------------
  list_quad1_poly<-lapply(1:length(gs),function(i){
    list_All_quads_poly[[i]]$quad.1.poly
  })
  list_quad2_poly<-lapply(1:length(gs),function(i){
    list_All_quads_poly[[i]]$quad.2.poly
  })
  list_quad3_poly<-lapply(1:length(gs),function(i){
    list_All_quads_poly[[i]]$quad.3.poly
  })
  list_quad4_poly<-lapply(1:length(gs),function(i){
    list_All_quads_poly[[i]]$quad.4.poly
  })
  list_quad1_pop<-lapply(1:length(gs),function(i){
    list_All_quads_poly[[i]]$quad.1
  })
  list_igm_poly<-lapply(1:length(gs),function(i){
    list_All_quads_poly[[i]]$igm.poly
  })
  list_quad.1<-lapply(1:length(gs),function(i){
    list_All_quads_poly[[i]]$quad.1
  })
  list_quad.2<-lapply(1:length(gs),function(i){
    list_All_quads_poly[[i]]$quad.2
  })
  list_quad.3<-lapply(1:length(gs),function(i){
    list_All_quads_poly[[i]]$quad.3
  })
  list_quad.4<-lapply(1:length(gs),function(i){
    list_All_quads_poly[[i]]$quad.4
  })
  names(list_quad1_poly)<-sampleNames(gs)
  node_IGigm<-add(gs,list_quad1_poly,parent="CD19+ B cells")
  names(list_quad2_poly)<-sampleNames(gs)
  node_IGigm<-add(gs,list_quad2_poly,parent="CD19+ B cells")
  names(list_quad3_poly)<-sampleNames(gs)
  node_IGigm<-add(gs,list_quad3_poly,parent="CD19+ B cells")
  names(list_quad4_poly)<-sampleNames(gs)
  node_IGigm<-add(gs,list_quad4_poly,parent="CD19+ B cells")
  names(list_igm_poly)<-sampleNames(gs)
  node_IGigm<-add(gs,list_igm_poly,parent="CD19+ B cells")
  recompute(gs)
  fs.quad1<-getData(gs,"IGD-IGM- B cells")
  fs.quad2<-getData(gs,"IGD+IGM- B cells")
  fs.quad3<-getData(gs,"IGD-IGM+ B cells")
  fs.quad4<-getData(gs,"IGD+IGM+ B cells")
  fs.igm<-getData(gs,"IGM")
  
  # show(autoplot(fs.19[[1]],"PE-CF594-A")),show(autoplot(fs.19[[1]],"FITC-A"))
  # show(autoplot(fs.19[[2]],"PE-CF594-A","FITC-A"))
  return(list(gs=gs,list_quad1_pop=list_quad1_pop,fs.igm=fs.igm,
              list_quad.1=list_quad.1,list_quad.2=list_quad.2,list_quad.3=list_quad.3,list_quad.4=list_quad.4))
  
}







# ------------------------Gating from CD20+ to CD10- -----------------------

gating_cd20_to_cd10neg<-function(fs.CD20pos,gs,fluorochrome.chans){
  print("gating_cd20_to_cd10neg")
  set.seed(40)
  list_CD10_poly<-lapply(1:length(gs), function(i){
    print(identifier(gs[[i]]@data[[i]]))
    
    cd10.gate <-deGate(fs.CD20pos[[i]],fluorochrome.chans["CD10"],upper=T,alpha=0.5)
    if(cd10.gate<1.8){
      group_string<-paste0("group_1:",identifier(gs[[1]]@data[[i]]))
      cd10.gate <-deGate(fs.CD20pos[[i]],fluorochrome.chans["CD10"],upper=T,use.upper=T,alpha=0.9)
    }else{
      group_string<-paste0("group_2:",identifier(gs[[1]]@data[[i]]))
    }
    print(paste0("CD10:",cd10.gate))
    cd10.neg<-flowDensity(fs.CD20pos[[i]],c(fluorochrome.chans["CD10"],fluorochrome.chans["CD27"]),position = c(F,NA),gates=c(cd10.gate,NA))
    cd10.pos<-flowDensity(fs.CD20pos[[i]],c(fluorochrome.chans["CD10"],fluorochrome.chans["CD27"]),position = c(T,NA),gates=c(cd10.gate,NA))
    cd10.pos_poly<-polygonGate(filterId = "CD10+",.gate=cd10.pos@filter)
    cd10.neg_poly<-polygonGate(filterId = "CD10-",.gate=cd10.neg@filter)
    
    return(list(cd10.pos_poly=cd10.pos_poly, cd10.neg_poly=cd10.neg_poly,groups=group_string,cd10.neg=cd10.neg,cd10.pos=cd10.pos))
  })
  #---------------------------
  list_CD10pos_poly<-lapply(1:length(gs),function(i){
    list_CD10_poly[[i]]$cd10.pos_poly
  })
  list_CD10neg_poly<-lapply(1:length(gs),function(i){
    list_CD10_poly[[i]]$cd10.neg_poly
  })
  list_groups_CD10<-lapply(1:length(gs),function(i){
    list_CD10_poly[[i]]$groups
  }) 
  list_cd10.pos<-lapply(1:length(gs),function(i){
    list_CD10_poly[[i]]$cd10.pos
  }) 
  list_cd10.neg<-lapply(1:length(gs),function(i){
    list_CD10_poly[[i]]$cd10.neg
  }) 
  names(list_CD10pos_poly)<-sampleNames(gs)
  names(list_CD10neg_poly)<-sampleNames(gs)
  node_CD10pos<-add(gs, list_CD10pos_poly, parent="CD19+CD20+")
  node_CD10neg<-add(gs, list_CD10neg_poly, parent="CD19+CD20+")
  recompute(gs)
  fs.cd10pos<-getData(gs,"CD10+")
  fs.cd10neg<-getData(gs,"CD10-")
  
  return(list(gs=gs,fs.cd10pos=fs.cd10pos,fs.cd10neg=fs.cd10neg,list_CD10pos_poly=list_CD10pos_poly,
              list_CD10neg_poly=list_CD10neg_poly,list_groups_CD10=list_groups_CD10,list_cd10.neg=list_cd10.neg,list_cd10.pos=list_cd10.pos))
}


# ------------------------ Gating from CD10- to Naive(CD27-,igd+) -----------------------
gating_cd10n_to_naive<-function(fs.cd10neg,gs,fluorochrome.chans){
  print("gating_cd10n_to_naive")
  set.seed(40)
  list_naive_poly<-lapply(1:length(gs),function(i){
    print(identifier(gs[[i]]@data[[i]]))
    # show(autoplot(fs.cd10neg[[1]],"PE-Cy7-A","PE-CF594-A"))
    # show(autoplot(fs.cd10neg[[1]],"PE-Cy7-A"))
    igd.gate <- deGate(fs.cd10neg[[i]],fluorochrome.chans["IgD"],upper=F,alpha=0.9)
    igd.gate<-igd.gate*0.95
    print(paste0("igd.gate:",igd.gate))
    #all_peaks <- fsApply(fs.10.neg,getPeaks,fluorochrome.chans["IgD"],tinypeak.removal=0.05)
    cd27.gate <- deGate(fs.cd10neg[[i]],fluorochrome.chans["CD27"])
    print(paste0("cd27.gate:",cd27.gate))
    if(cd27.gate<1.2){
      cd27.gate <- deGate(fs.cd10neg[[i]],fluorochrome.chans["CD27"],upper = T)
      cd27.gate_2 <- deGate(fs.cd10neg[[i]],fluorochrome.chans["CD27"],use.upper=T,upper = T,tinypeak.removal = 0.7,alpha = 0.001)
      if((cd27.gate_2-cd27.gate)>1){ # the possible gates are very far
        cd27.gate<-(cd27.gate+cd27.gate_2)/2# we take a solution in the middle
      }else{ # the possible gate are near each other
        cd27.gate<-cd27.gate_2
      }
    }
    print(paste0("cd27.gate:",cd27.gate))
    group_string<-paste0("group_1:",identifier(gs[[1]]@data[[i]]))
    quad.1 <- flowDensity(fs.cd10neg[[i]], c(fluorochrome.chans["CD27"],fluorochrome.chans["IgD"]),position=c(F,F),gates=c(cd27.gate*1.05,igd.gate))
    quad.2 <- flowDensity(fs.cd10neg[[i]],c(fluorochrome.chans["CD27"],fluorochrome.chans["IgD"]),position=c(T,F), gates=quad.1@gates)
    quad.3 <- flowDensity(fs.cd10neg[[i]], c(fluorochrome.chans["CD27"],fluorochrome.chans["IgD"]),position=c(F,T), gates=quad.1@gates)
    quad.4 <- flowDensity(fs.cd10neg[[i]], c(fluorochrome.chans["CD27"],fluorochrome.chans["IgD"]),position=c(T,T), gates=quad.1@gates)
    quad.1_poly<-polygonGate(filterId = "Atypical B cells",.gate=quad.1@filter)
    quad.2_poly<-polygonGate(filterId = "Unswitched memory B cells",.gate=quad.2@filter)
    quad.3_poly<-polygonGate(filterId = "Naive B cells",.gate=quad.3@filter)
    quad.4_poly<-polygonGate(filterId = "Switched memory B cells",.gate=quad.4@filter)
    
    return(list(quad.1_poly=quad.1_poly,quad.2_poly=quad.2_poly,quad.3_poly=quad.3_poly, quad.4_poly=quad.4_poly,
                quad.3=quad.3,groups=group_string,quad.1=quad.1,quad.2=quad.2,quad.4=quad.4))
  })
  #--------------------------------------------------
  list_quad1_poly<-lapply(1:length(gs),function(i){
    list_naive_poly[[i]]$quad.1_poly
  })
  list_quad2_poly<-lapply(1:length(gs),function(i){
    list_naive_poly[[i]]$quad.2_poly
  })
  list_quad3_poly<-lapply(1:length(gs),function(i){
    list_naive_poly[[i]]$quad.3_poly
  })
  list_quad4_poly<-lapply(1:length(gs),function(i){
    list_naive_poly[[i]]$quad.4_poly
  })
  list_quad.3<-lapply(1:length(gs),function(i){
    list_naive_poly[[i]]$quad.3
  })
  list_groups_Naive<-lapply(1:length(gs),function(i){
    list_naive_poly[[i]]$groups
  })
  list_quad.1<-lapply(1:length(gs),function(i){
    list_naive_poly[[i]]$quad.1
  })
  list_quad.2<-lapply(1:length(gs),function(i){
    list_naive_poly[[i]]$quad.2
  })
  list_quad.4<-lapply(1:length(gs),function(i){
    list_naive_poly[[i]]$quad.4
  })
  names(list_quad1_poly)<-sampleNames(gs)
  names(list_quad2_poly)<-sampleNames(gs)
  names(list_quad3_poly)<-sampleNames(gs)
  names(list_quad4_poly)<-sampleNames(gs)
  node_quad1<-add(gs, list_quad1_poly, parent="CD10-")
  node_quad2<-add(gs, list_quad2_poly, parent="CD10-")
  node_quad3<-add(gs, list_quad3_poly, parent="CD10-")
  node_quad4<-add(gs, list_quad4_poly, parent="CD10-")
  recompute(gs)
  fs.naive<-getData(gs,"Naive B cells")
  return(list(gs=gs,list_quad3_poly=list_quad3_poly,list_quad.3=list_quad.3,list_groups_Naive=list_groups_Naive,
              list_quad.1=list_quad.1,list_quad.2=list_quad.2,list_quad.4=list_quad.4))
}


# ------------------------ Gating from IGM+IGD+CD27- to Transitional cell-----------------------
gating_IGM_IGD_CD27_to_tran<-function(fs.igm,quad.3,gs,fluorochrome.chans){
  print("gating_IGM_IGD_CD27_to_tran")
  set.seed(40)
  list_transitional_poly<-lapply(1:length(gs),function(i){
    print(identifier(gs[[i]]@data[[i]]))
    cd27.igd <- flowDensity(fs.igm[[i]],  c(fluorochrome.chans["CD27"],fluorochrome.chans["IgD"]),position=c(F,T),gates=c(quad.3[[i]]@gates[1],quad.3[[i]]@gates[2]))
    gate.1 <- deGate(cd27.igd, fluorochrome.chans["CD10"],twin.factor = .9,tinypeak.removal =.9,upper=T, alpha=.5, magnitude=.1)
    gate.2 <- deGate(cd27.igd, fluorochrome.chans["CD38"],twin.factor = .9,tinypeak.removal =.9,upper=F, alpha=.9)
    print(paste0("gate.CD27:",gate.1))
    print(paste0("gate.CD38:",gate.2))
    if(gate.2<1){
      gate.2 <- deGate(cd27.igd, fluorochrome.chans["CD38"],twin.factor = .9,tinypeak.removal =.9,upper=T, alpha=0.1)
    }
    print(paste0("gate.CD38:",gate.2))
    group_string<-paste0("group_1:",identifier(gs[[1]]@data[[i]]))
    trans <- flowDensity(cd27.igd,c(fluorochrome.chans["CD10"],fluorochrome.chans["CD38"]),position=c(T,T),gates=c(gate.1,gate.2))
    trans.poly<-polygonGate(filterId = "Immature transition B cells",.gate=trans@filter)
    return(list(trans.poly=trans.poly,groups=group_string,trans=trans,f_cd27.igd=cd27.igd@flow.frame))
  })
  #------------------------------------
  list_trans_poly<-lapply(1:length(gs),function(i){
    list_transitional_poly[[i]]$trans.poly
  })
  list_groups_trans<-lapply(1:length(gs),function(i){
    list_transitional_poly[[i]]$groups
  })
  list_trans<-lapply(1:length(gs),function(i){
    list_transitional_poly[[i]]$trans
  })
  list_f_cd27.igd<-lapply(1:length(gs),function(i){
    list_transitional_poly[[i]]$f_cd27.igd
  })
  names(list_trans_poly)<-sampleNames(gs)
  node_trans<-add(gs, list_trans_poly, parent="IGM")
  recompute(gs)
  return(list(gs=gs,list_trans_poly=list_trans_poly,list_groups_trans=list_groups_trans,list_trans=list_trans,
              list_f_cd27.igd=list_f_cd27.igd))
} 

# ------------------------ Gating  from live to blast cells-----------------------
# show(autoplot(fs.live[[4]],"PE-A","SSC-A"))
# show(autoplot(fs.live,"SSC-A"))
# show(autoplot(fs.live[[4]],"PE-A"))
gating_live_to_blast<-function(fs.live,gs,fluorochrome.chans,scat.chans){
  print("gating_live_to_blast")
  set.seed(40)
  list_blast_poly<-lapply(1:length(gs),function(i){
    print(identifier(gs[[i]]@data[[i]]))
    cd.34_gate<- deGate(fs.live[[i]],fluorochrome.chans["CD34"],use.upper=T,upper=T,tinypeak.removal=0.5)
    print(paste0("cd.34_gate:",cd.34_gate))
    if(cd.34_gate>3){ # gate too much on right?
      all.peaks<-getPeaks(fs.live[[i]],fluorochrome.chans["CD34"])
      print(all.peaks)
      if(max(all.peaks$Peaks)>2.5){ # yes it is too much on the right, because there is a peak(probably the blast pop)
        # on the right and if it is >3 ,this gate does not catch this pop.
        cd.34_gate<-max(all.peaks$Peaks)*0.85
      }else if(max(all.peaks$Peaks)<2.5){ # No, the gate is correct, we only slighly put in on the left.
        cd.34_gate<-cd.34_gate*0.98
      }
    }else if(cd.34_gate<2){ # gate may be too much on left
      cd.34_gate<- deGate(fs.live[[i]],fluorochrome.chans["CD34"],use.upper=T,upper=T,alpha = 0.005)
      cd.34_gate<-cd.34_gate*1.15
    }
    ssca.gate <-deGate(fs.live[[i]],scat.chans["SSC-A"],bimodal = T)
    print(paste0("ssca.gate:",ssca.gate))
    print(paste0("cd.34_gate:",cd.34_gate))
    blast <- flowDensity(fs.live[[i]],c(fluorochrome.chans["CD34"],scat.chans["SSC-A"]),position = c(T,F),gates=c(cd.34_gate,ssca.gate))
    blast.poly<-polygonGate(filterId = "Blasts",.gate=blast@filter,blast=blast)
    return(list(blast.poly=blast.poly,blast=blast))
  })
  #--------------------------
  list_poly_blast<-lapply(1:length(gs),function(i){
    list_blast_poly[[i]]$blast.poly
  })
  list_groups_blast<-lapply(1:length(gs),function(i){
    list_blast_poly[[i]]$groups
  })
  list_blast<-lapply(1:length(gs),function(i){
    list_blast_poly[[i]]$blast
  })
  names(list_poly_blast)<-sampleNames(gs)
  node_blast<-add(gs, list_poly_blast, parent="Live cells")
  recompute(gs)
  return(list(gs=gs,list_poly_blast=list_poly_blast,list_groups_blast=list_groups_blast,list_blast=list_blast))
}
