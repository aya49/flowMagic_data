
#---------------- create margin gates -----------------------------
margin_gating <- function(fs,gs){
  set.seed(40)
  margin.gates <- fsApply(fs,removeMargins,c("FSC-A","SSC-A"),return.gate=T)
  rg <- lapply(1:nrow(margin.gates), function(x) return(rectangleGate(filterId="Margin", "FSC-A"=c(-1, margin.gates[x,1]),
                                                                      "SSC-A"=c(-1, margin.gates[x,2]))))
  names(rg) <-sampleNames(gs)
  nodeID1<-add(gs, rg)#,"Margin")
  recompute(gs)
  fs_margin<-getData(gs,"Margin")
  return(list(fs_margin=fs_margin,gs=gs))
}

#-------------- Compensation and Trasformation ----------------------
comp_and_transform <-function(gs,comp){
  set.seed(40)
  # ----- compensate and transform the margin population(all events population)
  gs<- compensate.flow(gs,comp = comp)
  recompute(gs)
  gs <-transform.flow(gs,remove.outliers=F,trans.chans = NULL)
  recompute(gs)
  fs.marg <- getData(gs, "Margin")
  return(list(fs.marg=fs.marg,gs=gs))
}    

#------ cleaning on margin events (all events) ---------
cleaning<-function(fs.marg,gs,n_cores,path.output){
  print("Cleaning")
  set.seed(40)  
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
  add(gs,clean.clust,parent = "/Margin", name = "Clean")
  recompute(gs)
  fs.clean <- getData(gs, "Clean_0")
  return(list(clean.inds=clean.inds,fs.clean=fs.clean,gs=gs))
}

# ------------------ Gating Time to get Singlets -------------------------------------------
gating_Time_to_singlets<-function(fs.clean,gs,pre_processing,n_cores,scat.chans){
  print("gating_Time_to_singlets")
  theta0 <- pi/5.2 # rotation angle
  set.seed(40)
  singlets_poly_list<-mclapply(1:length(gs),function(i){
    print(identifier(gs[[i]]@data[[i]]))
    rot <- rotate.data(fs.clean[[i]], c(scat.chans['FSC-A'], scat.chans['FSC-H']), theta = theta0)$data
    rot.gate <- deGate(rot, channel = scat.chans['FSC-H'], use.upper = T, upper = F, alpha = 0.05)
    print(paste0("rot.gate:",rot.gate))
    singlets <- rotate.fd(flowDensity(rot, channels = c(scat.chans['FSC-A'], scat.chans['FSC-H']), position = c(NA, T), gates = c(NA, rot.gate)),angle = theta0)@filter
    singlets_pop <- rotate.fd(flowDensity(rot, channels = c(scat.chans['FSC-A'], scat.chans['FSC-H']), position = c(NA, T), gates = c(NA, rot.gate)),angle = theta0)
    
    sngl.poly <- polygonGate(filterId = "Singlets",.gate=singlets)
    return(list(sngl.poly=sngl.poly,singlets=singlets,singlets_pop=singlets_pop))
  },mc.cores=n_cores,mc.set.seed = F)
  list_singlets_poly <- lapply(1:length(gs),function(i){
    singlets_poly_list[[i]]$sngl.poly
  })
  list_singlets <- lapply(1:length(gs),function(i){
    singlets_poly_list[[i]]$singlets
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
  # estraiamo popolazione dei Singlets
  fs.sngl<-getData(gs,"Singlets")
  
  return(list(fs.sngl=fs.sngl,gs=gs,list_singlets=list_singlets,list_singlets_pop=list_singlets_pop))
}


# --------------------------------- Gating singlets to get Beads, Size -------------------------------------
gating_singlets_to_Size_beads<-function(fs.sngl,gs,n_cores,scat.chans){
  print("gating_Size_beads")
  # show(autoplot(fs.sngl[[1]],"FSC-A","SSC-A"))
  set.seed(40)
  beads_size_poly_list<-mclapply(1:length(gs),function(i){
    print(identifier(gs[[i]]@data[[i]]))
    theta0 <- (1/12)*pi # defining rotation angle (15 degree rotation)
    rot <- rotate.data(fs.sngl[[i]], c(scat.chans['FSC-A'], scat.chans['SSC-A']), theta = theta0)$data
    #show(autoplot(rot,"FSC-A","SSC-A"))
    # show(autoplot(rot,"SSC-A"))
    all.cuts<-deGate(rot, channel = c(scat.chans['SSC-A']), all.cuts = T,tinypeak.removal = 0.01)
    all.peaks<-getPeaks(rot, channel = c(scat.chans['SSC-A']), tinypeak.removal = 0.01)
    print(all.peaks)
    print(paste0("all.cuts:",all.cuts))
    inds<-which(all.cuts<145000)
    all.cuts<-all.cuts[inds]
    max_cuts_gate<-max(all.cuts)
    inds<-which(all.peaks$Peaks<110000)
    highest_not_bead_peak<-max(all.peaks$Peaks[inds]) # it is the peak of the size gate,just below the beads gate
    if(max_cuts_gate<=50000 && max(all.peaks$Peaks)>190000){ # beads peak is very high and cuts gate is too much low
      all.cuts<-deGate(rot, channel = c(scat.chans['SSC-A']), all.cuts = T,tinypeak.removal = 0.01)
      inds<-which(all.cuts<170000)
      all.cuts<-all.cuts[inds]
      max_cuts_gate<-max(all.cuts)
    }
    print(paste0("max_cuts_gate:",max_cuts_gate))
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
    if(highest_not_bead_peak>70000){ # size peak very near the beads population
      print("WARNING: rotation angle must change! size pop is too near bead pop")
      theta0 <- (1/6)*pi # defining rotation angle (30 degree rotation)
      rot <- rotate.data(fs.sngl[[i]], c(scat.chans['FSC-A'], scat.chans['SSC-A']), theta = theta0)$data
      all.peaks<-getPeaks(rot, channel = c(scat.chans['SSC-A']), tinypeak.removal = 0.01)
      print(all.peaks)
      #show(autoplot(rot,"FSC-A","SSC-A"))
      inds<-which(all.peaks$Peaks>80000)
      bead_peak_rotated<-max(all.peaks$Peaks[inds]) # it is the peak of the bead gate
      max_cuts_gate<-bead_peak_rotated*0.75
      if(highest_not_bead_peak<80000){
        max_cuts_gate<-bead_peak_rotated*0.68
      }
      print(paste0("max_cuts_gate_rotated:",max_cuts_gate))
    }
    temp<- flowDensity(rot, c(scat.chans['FSC-A'], scat.chans['SSC-A']), position = c(NA, T), gates = c(NA,max_cuts_gate))
    #show(autoplot(getflowFrame(temp),"FSC-A","SSC-A"))
    #show(autoplot(getflowFrame(temp),"SSC-A"))
    rot.SSCA.gate<-deGate(temp, channel = c(scat.chans['SSC-A']), use.upper=T,upper=F,tinypeak.removal = 0.01,alpha = 0.001)
    if(max_cuts_gate>=135000){ 
      rot.SSCA.gate<-rot.SSCA.gate*1.05
    }
    if(rot.SSCA.gate<40000){
      rot.SSCA.gate<-deGate(temp, channel = c(scat.chans['SSC-A']), use.upper=T,upper=F,tinypeak.removal = 0.01,alpha = 0.3)
      
    }
    print(paste0("rot.SSCA.gate:",rot.SSCA.gate))
    beads<- rotate.fd(flowDensity(rot, c(scat.chans['FSC-A'], scat.chans['SSC-A']), position = c(NA, T), gates = c(NA,rot.SSCA.gate)),angle = theta0)
    # show(autoplot(getflowFrame(beads),"FSC-A","SSC-A"))
    # show(autoplot(getflowFrame(beads),"FSC-A"))
    #beads <- flowDensity(rot, channels = c(scat.chans['FSC-A'], scat.chans['SSC-A']), position = c(NA, T), gates = c(NA, rot.SSCA.gate*1.05))
    FSCAbead.gate<-deGate(beads, channel = c(scat.chans['SSC-A']), upper=T,alpha = 0.8)
    print(paste0("FSCAbead.gate:",FSCAbead.gate))
    if((FSCAbead.gate<70000)||(FSCAbead.gate>100000)){
      all.peaks<-getPeaks(beads, channel = c(scat.chans['FSC-A']),tinypeak.removal = 0.1)
      inds<-which(all.peaks$Peaks>40000)
      all.peaks$Peaks<-all.peaks$Peaks[inds]
      FSCAbead.gate<-min(all.peaks$Peaks)
      FSCAbead.gate<-FSCAbead.gate*1.15
      print(FSCAbead.gate)
    }
    
    print(paste0("FSCAbead.gate:",FSCAbead.gate))
    beads.poly<- polygonGate(filterId = "Beads",.gate=flowDensity(beads, channels = c(scat.chans['FSC-A'], scat.chans['SSC-A']), 
                                                                  position = c(F, NA), gates = c(FSCAbead.gate, NA), ellip.gate = T, scale = 0.96)@filter)
    beads<-flowDensity(beads, channels = c(scat.chans['FSC-A'], scat.chans['SSC-A']),position = c(F, NA), gates = c(FSCAbead.gate, NA), ellip.gate = T, scale = 0.96)
    # ------- find the size population 

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
    if(left_peak>140000){ # left peak on the extreme right (maybe too many low cells in the file)
      rot.gate <- deGate(temp, channel = scat.chans['FSC-A'],use.upper=F, upper = F)
      if(rot.gate>100000){# gate still incorrect
        rot.gate<-left_peak/2*0.90
      }
      print(paste0("temp_left_peak_100000:",rot.gate))
    }else if((left_peak>=47000)&&(left_peak<128000)){ # it means that deGate and so getPeaks don't see the first peak or it is absent
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
      rot.gate<-deGate(rot, channel = scat.chans['FSC-A'],after.peak = T,alpha = 0.9)
      if(rot.gate>65000){ # gate too much on the right
        rot.gate<-(rot.gate+left_peak)/2
      }
    }
    
    
    print(paste0("rot.gate:",rot.gate))
    if((rot.gate>85000)&&(left_peak<18000)){
      rot.gate<-left_peak*4
    }
    
    print(paste0("rot.gate:",rot.gate))
    size.temp<- rotate.fd(flowDensity(rot, c(scat.chans['FSC-A'], scat.chans['SSC-A']), position = c(T, NA), gates = c(rot.gate, NA)),angle = theta0)
    size.poly <- polygonGate(filterId = "All cells",.gate=size.temp@filter)
    
    return(list(beads.poly=beads.poly,size.poly=size.poly,beads=beads,size=size.temp))
  },mc.cores = n_cores, mc.set.seed = F)
  beads_poly_list<-lapply(1:length(gs),function(i){
    beads_size_poly_list[[i]]$beads.poly
  })
  size_poly_list<-lapply(1:length(gs),function(i){
    beads_size_poly_list[[i]]$size.poly
  })
  
  list_beads<-lapply(1:length(gs),function(i){
    beads_size_poly_list[[i]]$beads
  })
  list_size<-lapply(1:length(gs),function(i){
    beads_size_poly_list[[i]]$size
  })
  names(beads_poly_list)<-sampleNames(gs)
  nodeID0<-add(gs,beads_poly_list,parent="Singlets")
  recompute(gs)
  names(size_poly_list)<-sampleNames(gs)
  names(list_size)<-sampleNames(gs)
  nodeID0<-add(gs,size_poly_list,parent="Singlets")
  recompute(gs)
  fs.size<-getData(gs,"All cells")
  return(list(fs.size=fs.size,gs=gs,beads_poly_list=beads_poly_list,
              size_poly_list=size_poly_list,list_beads=list_beads,list_size=list_size))
}

# ------------------ gating to live from the Size population (All cells population) ---------
gating_size_to_live<-function(fs.size,gs,n_cores,fluorochrome.chans,scat.chans){
  print("gating_live")
  set.seed(40)
  live_poly_list<-mclapply(1:length(gs),function(i){
    print(identifier(gs[[1]]@data[[i]]))
    # in order to obtain an oblique gates,we need to rotate the data,
    # but we cannot rotate using big angles,because the data are scaled 
    # in a different way along the x axis and y axis,
    # I use the min.max option of the rotate.data function 
    #to find the right angle rotation, and 
    # I play around this angle(only a little bit,the numbers must be very small)
    output <- rotate.data(fs.size[[i]], c(fluorochrome.chans["Live"], scat.chans['SSC-A']), min.max = T)
    theta0 <- output$theta
    theta0<-theta0/8
    rot <- rotate.data(fs.size[[i]], c(fluorochrome.chans["Live"], scat.chans['SSC-A']), theta = theta0)$data
    #show(autoplot(rot,"APC-eF780-A","SSC-A"))
    #show(autoplot(rot,"APC-eF780-A"))
    rot.live.gate_all.cuts <- deGate(rot,"APC-eF780-A",all.cuts=T)
    all.peaks <- getPeaks(rot,"APC-eF780-A",tinypeak.removal = 0.1)
    rot.live.gate_all.peaks<-all.peaks$Peaks
    print(paste0("all.peaks.rot:",rot.live.gate_all.peaks))
    print(paste0("all.cuts.rot:",rot.live.gate_all.cuts))
    max_peak_rot<-max(rot.live.gate_all.peaks)
    print(paste0("max_peak_rot:",max_peak_rot))
    if(max_peak_rot>=2.82){ # the max peak is on the right,we put the gate on the left of this peak
      rot.temp.live.gate<-max_peak_rot
      temp<-flowDensity(rot, c(fluorochrome.chans['Live'], scat.chans['SSC-A']), position = c(F, NA), gates = c(max_peak_rot, NA))
      rot.live.gate<-deGate(temp, channel = fluorochrome.chans['Live'], use.upper=T, upper = T, alpha = 0.04)
      #show(autoplot(getflowFrame(temp),"APC-eF780-A","SSC-A"))
      #show(autoplot(getflowFrame(temp),"APC-eF780-A"))
      print(rot.live.gate)
      if(rot.live.gate>=2.9){
        rot.live.gate<-deGate(temp, channel = fluorochrome.chans['Live'], use.upper=T,upper = T, alpha = 0.01,tinypeak.removal = 0.5)
      }
      
      if((rot.live.gate>=3.2)||(rot.live.gate<2)){
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
      }
      
      
    }else{ # the max peak is more on the left,we put the gate on the right of this peak
      rot.live.gate<-deGate(rot, channel = fluorochrome.chans['Live'], use.upper=T, upper = T, alpha = 0.01)
      print(paste0("rot.live.gate_0.01:",rot.live.gate))
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
        inds<-which(rot.live.gate_all.cuts<3.0)
        rot.live.gate_all.cuts<-rot.live.gate_all.cuts[inds]
        rot.live.gate<-max(rot.live.gate_all.cuts)
        print(paste0("rot.live.gate_cuts:",rot.live.gate))
      }
      #-------------------------------- 
      # gate too much on the left
      if((rot.live.gate<1.7)&&(max_peak_rot>=1)){
        rot.live.gate<-deGate(rot, channel = fluorochrome.chans['Live'], use.upper=T, upper = T, alpha = 0.2,tinypeak.removal = 0.3)
        print(paste0("rot.live.gate_0.2:",rot.live.gate))
      } else if((max_peak_rot<1)&&(rot.live.gate<1.7)){
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
    return(list(live.poly=live.poly,live=live))
  },mc.cores = n_cores,mc.set.seed = F)
  
  list_live_poly<-lapply(1:length(gs),function(i){
    live_poly_list[[i]]$live.poly
  })
  list_live<-lapply(1:length(gs),function(i){
    live_poly_list[[i]]$live
  })
  names(list_live_poly)<-sampleNames(gs)
  names(list_live)<-sampleNames(gs)
  nodeID0<-add(gs,list_live_poly,parent="All cells")
  recompute(gs)
  fs.live <- getData(gs, "Live cells")
  return(list(fs.live=fs.live,gs=gs,live_poly_list=list_live_poly,list_live=list_live))
}

# ---------------------- Gating Live to get Granulocytes e CD45+CD66- -------------------------------------------
gating_live_to_granulocytes<-function(fs.live,gs,n_cores,fluorochrome.chans){
  print("gating_granulocytes")
  # show(autoplot(fs.live[[1]],"BV711-A","V450-A"))
  set.seed(40)
  list_poly_gran_cd45<-mclapply(1:length(gs),function(i){
    # CD66.gate <- deGate(fs.live[[i]], channel = c(fluorochrome.chans['CD66']), all.cuts = T, adjust.dens = .8)[1]
    CD66.gate <- deGate(fs.live[[i]], channel = c(fluorochrome.chans['CD66']), upper = T)[1]
    print(paste0("CD66.gate:",CD66.gate))
    temp <- flowDensity(fs.live[[i]] , channels = c(fluorochrome.chans['CD66'], fluorochrome.chans['CD45']), position = c(T, NA), gates = c(CD66.gate, NA))
    # plotDens(fs.live[[1]],c(fluorochrome.chans['CD66'],fluorochrome.chans['CD45']))
    CD66.gate.up<- deGate(temp,fluorochrome.chans['CD66'],upper=F)
    if(CD66.gate.up>3){
      CD66.gate.up <- deGate(temp,fluorochrome.chans['CD66'],upper=F,use.upper=T)
    }
    if(CD66.gate.up<CD66.gate){
      CD66.gate.up<-CD66.gate*1.02
    }
    if(CD66.gate<=1){
      CD66.gate<-CD66.gate.up*0.90
    }
    print(paste0("CD66.gate.up:",CD66.gate.up))
    granulocytes<- flowDensity(fs.live[[i]], channels = c(fluorochrome.chans['CD66'], fluorochrome.chans['CD45']), position = c(T, NA), gates = c(CD66.gate.up, NA))
    CD45granulo.upper.gate <- deGate(granulocytes, channel = c(fluorochrome.chans['CD45']), use.upper = T, upper = T, alpha = 0.05) + 0.05
    CD45granulo.lower.gate <- deGate(granulocytes, channel = c(fluorochrome.chans['CD45']), use.upper = T, upper = F,alpha=0.01)
    granulocytes <- flowDensity(granulocytes, channels = c(fluorochrome.chans['CD66'], fluorochrome.chans['CD45']), position = c(NA, T), gates = c(NA, CD45granulo.lower.gate))
    granulocytes <- flowDensity(granulocytes, channels = c(fluorochrome.chans['CD66'], fluorochrome.chans['CD45']), position = c(NA, F), gates = c(NA, CD45granulo.upper.gate))

    #------- now we are finding the CD45+CD66- pop
    CD66negCD45pos<-flowDensity(fs.live[[i]], channels = c(fluorochrome.chans['CD66'], fluorochrome.chans['CD45']), position = c(F, NA), gates = c(CD66.gate, NA))
    CD45.gate <- deGate(CD66negCD45pos, channel = c(fluorochrome.chans['CD45']), upper = F, alpha = 0.3)
    print(paste0("CD45.gate:",CD45.gate))
    if(CD45.gate<2){
      CD45.gate <- deGate(CD66negCD45pos, channel = c(fluorochrome.chans['CD45']), after.peak = T,alpha = 0.3)
    }
    print(paste0("CD45.gate:",CD45.gate))
    
    CD66negCD45pos<- flowDensity(fs.live[[i]], channels = c(fluorochrome.chans['CD66'], fluorochrome.chans['CD45']), position = c(F, T), gates = c(CD66.gate, CD45.gate))
    gran.poly <- polygonGate(filterId = "Granulocytes",.gate=granulocytes@filter)
    cd66.45.poly <- polygonGate(filterId = "CD45+CD66- (Non Granulocytes)",.gate  =CD66negCD45pos@filter)
    return(list(gran.poly=gran.poly,cd66.45.poly=cd66.45.poly,CD66negCD45pos=CD66negCD45pos,granulocytes=granulocytes))
  },mc.cores=n_cores,mc.set.seed = F)
  
  lista_poly_gran<-lapply(1:length(gs),function(i){
    list_poly_gran_cd45[[i]]$gran.poly
  })
  lista_poly_cd45<-lapply(1:length(gs),function(i){
    list_poly_gran_cd45[[i]]$cd66.45.poly
  })
  lista_CD66negCD45pos<-lapply(1:length(gs),function(i){
    list_poly_gran_cd45[[i]]$CD66negCD45pos
  })
  lista_granulocytes<-lapply(1:length(gs),function(i){
    list_poly_gran_cd45[[i]]$granulocytes
  })
  names(lista_poly_gran)<-sampleNames(gs)
  names(lista_CD66negCD45pos)<-sampleNames(gs)
  names(lista_granulocytes)<-sampleNames(gs)
  
  nodeID0<-add(gs,lista_poly_gran,parent="Live cells")
  names(lista_poly_cd45)<-sampleNames(gs)
  nodeID0<-add(gs,lista_poly_cd45,parent="Live cells")
  recompute(gs) 
  # Estraiamo il flowSet dei granulociti.
  fs.gran <- getData(gs, "Granulocytes")
  fs.66n45p<-getData(gs,"CD45+CD66- (Non Granulocytes)")
  return(list(fs.gran=fs.gran,gs=gs,fs.66n45p=fs.66n45p,lista_poly_gran=lista_poly_gran,lista_poly_cd45=lista_poly_cd45,lista_CD66negCD45pos=lista_CD66negCD45pos,lista_granulocytes=lista_granulocytes))
}


# -----------------------------  Gating CD45+CD66- to obtain monocytes(HLADR+ CD14+) e HLADR+ CD14-  ---------------------------------------------------------
gating_ngran_to_DRpos<-function(fs.66n45p,gs,fluorochrome.chans){
  print("Gating_DRpos")
  set.seed(40)
  lista_poly_mono_drp<-lapply(1:length(gs),function(i){
    print(identifier(gs[[i]]@data[[i]]))
    CD14.gate <- deGate(fs.66n45p[[i]], channel = c(fluorochrome.chans['CD14']), upper=T, bimodal = T, twin.factor = 0.8,alpha=.2,magnitude = .7,tinypeak.removal = 0.01)
    #print(CD14.gate)
    all.peaks<-getPeaks(fs.66n45p[[i]], channel = c(fluorochrome.chans['CD14']),tinypeak.removal = 0.15)
    print(all.peaks)
    ind<-which(all.peaks$P.h==max(all.peaks$P.h))
    intense_peak<-all.peaks$Peaks[ind]
    #-----------------------------------------------
    if(length(all.peaks$Peaks)==1){
      if(intense_peak>=2.5){
        if(intense_peak<=2.65){
          CD14.gate <- intense_peak
        }else{
          if(intense_peak<2.80){
            CD14.gate <- intense_peak*0.97
          }else if(intense_peak>=2.80){
            CD14.gate <- intense_peak*0.90
          }
        }
        
      }else if(intense_peak<2.5){
        if(intense_peak>2.2){
          CD14.gate <- intense_peak*1.07
        }else if(intense_peak>=1.9){
          CD14.gate <- intense_peak*1.15
        }else if(intense_peak>=1.5){
          CD14.gate <- intense_peak*1.30
        }else{
          CD14.gate <- intense_peak*1.40
        }
      }
    }
    #-----------------------------------------------
    if(length(all.peaks$Peaks)==2){
      diff<-all.peaks$Peaks[2]-all.peaks$Peaks[1]
      print(paste0("diff:",diff))
      if(diff<=0.3){
        # we consider them as one peak. 
        if(intense_peak<=2){
          #Below 2 so we gate in the same way as the second first peak situation
          # but we cannot use the same cose but the second almost twin peak alters the result
          CD14.gate <- all.peaks$Peaks[2]*1.08
        }else if(intense_peak>2){
          # above 2.1 same situation as the first for loop
          CD14.gate <- intense_peak*0.80
        }
      }else if(diff>0.3){
        # we consider them as two different peaks.
        #----------- the second peak is in a very  high position
        if(all.peaks$Peaks[2]>2.5){
          if(all.peaks$Peaks[1]>=2.15){ # the first peak is near the second one(both peaks in a very high position)
            CD14.gate<-all.peaks$Peaks[1] # gate in the first peak
            if(diff>=0.4){ #  the first peak is near the second one but not so near. Gate in the middle
              CD14.gate<-(all.peaks$Peaks[1]+all.peaks$Peaks[2])/2  
            }
            if(diff>1){ # the two peaks are at the limit between being near and being far
              CD14.gate<-all.peaks$Peaks[1]*1.10 # gate slightly upper the first peak
            }
            if((all.peaks$Peaks[1]>=2.3)&&(diff<0.4)){ # the first peak is very near the second one. Despite we consider them as two different peaks. We put the gate slightly under the second peak
              CD14.gate<-all.peaks$Peaks[1]*0.95 # sligtly under the first peak
            }
          }else if(all.peaks$Peaks[1]<2.15){ # the first peak is far from the second peak
            CD14.gate<-deGate(fs.66n45p[[i]], channel = c(fluorochrome.chans['CD14']), upper=T, tinypeak.removal = 0.15) # gate in the middle of the peaks
            if(all.peaks$Peaks[1]<0.9){ # the first peak is too far from the second peak to put the gate in the middle
              CD14.gate<-all.peaks$Peaks[2]*0.90 # gate lower than the second peak
            }
          }
          #-------- the second peak is in a lower position (but not too low)
        }else if(all.peaks$Peaks[2]<2.5 && all.peaks$Peaks[2]>2){ 
          CD14.gate<-deGate(fs.66n45p[[i]], channel = c(fluorochrome.chans['CD14']), upper=T, tinypeak.removal = 0.15) # gate in the middle of the peaks
          if(all.peaks$Peaks[1]<0.9){ # the first peak is too far from the second peak to put the gate in the middle
            CD14.gate<-all.peaks$Peaks[2]*0.90 # gate lower than the second peak
          }
          #-------- the second peak is in a very low position
        }else if(all.peaks$Peaks[2]<2){
          CD14.gate<-deGate(fs.66n45p[[i]], channel = c(fluorochrome.chans['CD14']), use.upper=T, upper=T,tinypeak.removal = 0.15) # gate after second peak
          diff_2<-CD14.gate-all.peaks$Peaks[2]
          if(diff_2>0.3){
            CD14.gate<-all.peaks$Peaks[2]*1.20
          }
        }
      }
    }
    #--------------------------------
    if(length(all.peaks$Peaks)==3){
      if(CD14.gate>1.8){
        CD14.gate<-CD14.gate
        diff<-all.peaks$Peaks[2]-all.peaks$Peaks[1]
        if((all.peaks$Peaks[3])>=2.8 && diff>0.3 && (all.peaks$Peaks[2])>=2 && (all.peaks$Peaks[1]>=2)){ # third peak far away,the second and first peak in high position and distant
          CD14.gate<-(all.peaks$Peaks[2]+all.peaks$Peaks[1])/2 # I put the gate between the second and first peak
        }
      }else{
        if((all.peaks$Peaks[2])>=2 && (all.peaks$Peaks[1]>=2)){
          CD14.gate<-(all.peaks$Peaks[2]+all.peaks$Peaks[1])/2 # I put the gate between the second and first peak
        }else if((all.peaks$Peaks[2])>=2 && (all.peaks$Peaks[1]<=2)){
          if((all.peaks$Peaks[2])>2.3){
            CD14.gate<-all.peaks$Peaks[2]*0.90 # under the second peak
          }else{
            CD14.gate<-all.peaks$Peaks[2]*1.10 # above the second peak
          }
        }
        if(CD14.gate<1){
          inds<-which(all.peaks$Peaks>1)
          all.peaks$Peaks<-all.peaks$Peaks[inds]
          n<-length(all.peaks$Peaks)
          CD14.gate<-sum(all.peaks$Peaks)/n
        }
      }
    }
    
    print(paste0("CD14.gate:",CD14.gate))
    CD14neg <- flowDensity(fs.66n45p[[i]], channels = c(fluorochrome.chans['HLA-DR'], fluorochrome.chans['CD14']), position = c(NA, F), gates = c(NA, CD14.gate))
    #show(autoplot(CD14neg@flow.frame,"eF605-A","V500-A"))
    HLADR.gate <- c(deGate(CD14neg, channel = c(fluorochrome.chans['HLA-DR']), all.cuts = T, upper = T, percentile = NA),
                    deGate(CD14neg, channel = c(fluorochrome.chans['HLA-DR']), use.upper = T, upper = T))
    maxDens <- density(getflowFrame(CD14neg)@exprs[,c(fluorochrome.chans['HLA-DR'])])
    start.idx <- which.min(abs(maxDens$x - 1.5))
    HLADR.gate <- min(HLADR.gate[which(HLADR.gate > maxDens$x[which.max(maxDens$y[start.idx:length(maxDens$y)]) + start.idx - 1])])
    all.peaks<-getPeaks(CD14neg, channel = fluorochrome.chans['HLA-DR'], tinypeak.removal = 0.4)
    print(all.peaks)
    inds<-which(all.peaks$Peaks<3.1)
    last_HLADR_peak<-max(all.peaks$Peaks[inds])
    if(last_HLADR_peak=="-Inf"){
      last_HLADR_peak<-max(all.peaks$Peaks)
    }
    first_HLADR_peak<-min(all.peaks$Peaks[inds])
    print(paste0("last_HLADR_peak:",last_HLADR_peak))
    print(paste0("HLADR.gate:",HLADR.gate))
    #show(autoplot(fs.66n45p[[1]],"eF605-A","V500-A"))
    #print(HLADR.gate)
    if((HLADR.gate<=2)||(HLADR.gate>3.1)){
      if(last_HLADR_peak<2.2){ # the max hladr peak is on the left,the gate should be on extremely on the right
        HLADR.gate<-deGate(CD14neg, channel = c(fluorochrome.chans['HLA-DR']), use.upper = T, upper = T,alpha=0.01,tinypeak.removal = 0.4)
        if(HLADR.gate>=3.2){
          HLADR.gate<-last_HLADR_peak*1.20
          if(last_HLADR_peak<2.0){
            HLADR.gate<-last_HLADR_peak*1.30
          }
        }
        #print(HLADR.gate)
      }else if(last_HLADR_peak<2.9){ # the max hladr peak is on the rigth but only sligtly, the gate should be still on the right, despite not too much
        HLADR.gate<-deGate(CD14neg, channel = c(fluorochrome.chans['HLA-DR']), use.upper = T, upper = T, alpha=0.5, tinypeak.removal = 0.4)
      }else{ 
        all.cuts<-deGate(CD14neg, channel = c(fluorochrome.chans['HLA-DR']), all.cuts = T)
        print(paste0("all.cuts:",all.cuts))
        inds<-which(all.cuts<=last_HLADR_peak)
        if(length(inds)!=0){
          all.cuts<-all.cuts[inds]
          HLADR.gate<-max(all.cuts)
        }
      }
    }
    print(paste0("HLADR.gate:",HLADR.gate))
    print(abs(HLADR.gate-last_HLADR_peak))
    if(abs(HLADR.gate-last_HLADR_peak)<0.3){ # the gate is too near the peak
      if(last_HLADR_peak<HLADR.gate){
        HLADR.gate<-HLADR.gate*1.10
        print(paste0("HLADR.gate:",HLADR.gate))
      }else{
        HLADR.gate<-HLADR.gate*0.90
        print(paste0("HLADR.gate:",HLADR.gate))
      }
    }
    print(paste0("HLADR.gate:",HLADR.gate))
    DRn.14n <-flowDensity(CD14neg, channels = c(fluorochrome.chans['HLA-DR'], fluorochrome.chans['CD14']), position = c(F, NA), gates = c(HLADR.gate, NA))

    # ---- get the monocites pop from the CD66-CD45+
    # notSubFrame = Remove a subset of a FlowFrame object specified by gates from the flowDensity method.
    temp <- notSubFrame(fs.66n45p[[i]], channels = c(fluorochrome.chans['HLA-DR'], fluorochrome.chans['CD14']), position = c(F, NA),filter= DRn.14n@filter)

    rot <- rotate.data(temp@flow.frame,c(fluorochrome.chans['HLA-DR'], fluorochrome.chans['CD14']),theta = -pi/4)

    #show(ggplot(rot$data, aes(x = `eF605-A`, y = `V500-A`)) + geom_hex(bins = 128) + scale_fill_gradientn(colours = gray.colors(9)))

    all.peaks<-getPeaks(rot$data, channel = fluorochrome.chans['HLA-DR'], tinypeak.removal = 0.3)
    print(all.peaks)
    if(length(all.peaks$Peaks)>1){
      ########################## multiple peaks ##################################33
      #print(hladr.gate)
      last_peak<-max(all.peaks$Peaks)
      print(paste0("last_peak:",last_peak))
      ###########################################################
      # the last peak is on the extreme right
      if(last_peak>0.2){
        # we need to position the gate on the left of last peak
        min_temp_peak<-min(all.peaks$Peaks)
        max_temp_peak<-max(all.peaks$Peaks)
        hladr.gate <- deGate(rot$data,fluorochrome.chans['HLA-DR'],upper = F,alpha=0.5,tinypeak.removal = 0.3)
        print(hladr.gate)
        if(last_HLADR_peak>2.5 && last_HLADR_peak<2.7){
          hladr.gate<-hladr.gate*1.22
          if(hladr.gate<0.5){
            hladr.gate<-(hladr.gate+last_peak)/2
          }
        }
        print(paste0("hladr.gate_1:",hladr.gate))
        #---------------------------------------------------------
        # the gate is too much on the left
        if(hladr.gate<(-0.25)){ 
          hladr.gate <- deGate(rot$data,fluorochrome.chans['HLA-DR'], use.upper = T, upper = T,alpha=0.05,tinypeak.removal = 0.3)
          hladr.gate_1 <- deGate(rot$data,fluorochrome.chans['HLA-DR'], use.upper = T, upper = T,alpha=0.05)
          hladr.gate<-max(hladr.gate_1,hladr.gate)
          print(paste0("hladr.gate_01:",hladr.gate))
          print(hladr.gate-min_temp_peak)
          if(((hladr.gate-min_temp_peak)>0.5)&&(max_temp_peak>0.2)){ # with max peak on the right,the gate is too much distant from the min peak. The gate should be between the min peak and the max peak of the temp population
            hladr.gate<- (min_temp_peak+max_temp_peak)/2
            if(min_temp_peak<(-1)){ # but if min peak is too much on the left,the mean is biased too much on the left,so we calculated the gate in another way
              hladr.gate<-max_temp_peak/2
            }
            if(hladr.gate<(-0.25)){ # if the gate is still not correct
              hladr.gate<-(last_peak+hladr.gate)/2
            }
          }else if(((hladr.gate-max_temp_peak)>0.5)&&(max_temp_peak<0.2)){ # if max peak is on the left,and the gate is too much distant from this peak
            hladr.gate <- deGate(rot$data,fluorochrome.chans['HLA-DR'], use.upper = T, upper = T,alpha=0.9,tinypeak.removal = 0.3)
            if((hladr.gate-max_temp_peak)>0.5){
              hladr.gate<-(hladr.gate+max_temp_peak)/2
            }
            if(hladr.gate<(-0.25)){ # if the gate is still not correct
              hladr.gate<-(last_peak+hladr.gate)/2
            }
          }
        }
        print(paste0("hladr.gate_2:",hladr.gate))
        #--------------------------------------------------------------------
        # the gate is too much on the right
        if((hladr.gate>0.6)&&(last_HLADR_peak<2.5)&&(last_peak<0.3)){ 
          if(length(all.peaks$Peaks)==2){
            hladr.gate <- deGate(rot$data,fluorochrome.chans['HLA-DR'], upper = T,alpha=0.9, tinypeak.removal = 0.5) # we test the gate without removing peaks using upper=T to put the gate between the peaks.
            print(paste0("temp:",hladr.gate))
            if(hladr.gate>0.3){ # if the gate is still too much on the right
              if(all.peaks$Peaks[1]>(-0.8)){
                hladr.gate <- (all.peaks$Peaks[1]+all.peaks$Peaks[2])/2 # we test the gate without removing peaks using upper=T to put the gate between the peaks.
              }else{
                hladr.gate <- deGate(rot$data,fluorochrome.chans['HLA-DR'],use.percentile = T,percentile = 0.9) # higher the percentlile,more on the right the gate will be
              }
            }else if(hladr.gate<(-0.1)){ # gate again too much on the left
              hladr.gate<-max(all.peaks$Peaks)*0.80
            }
          }else{
            hladr.gate <- deGate(rot$data,fluorochrome.chans['HLA-DR'],upper = T,twin.factor = 0.8)
            hladr.gate_1 <- deGate(rot$data,fluorochrome.chans['HLA-DR'],upper = T,bimodal = T) 
            hladr.gate<-max(hladr.gate_1,hladr.gate)
            print(paste0("temp:",hladr.gate))
            if((hladr.gate>0.3)&&(max(all.peaks$Peaks)>0.1)){ # gate is still too much on the right,it should be on the left of the peak
              hladr.gate <- deGate(rot$data,fluorochrome.chans['HLA-DR'],upper = F,alpha = 0.9)
              if(hladr.gate<(-0.1)){ # gate again too much on the left
                hladr.gate<-max(all.peaks$Peaks)*0.80
              }
            }else if((hladr.gate>0.4)&&(max(all.peaks$Peaks)<0.1)){ # gate is still too much on the right,but must be still on the right of peak
              hladr.gate <- (all.peaks$Peaks+hladr.gate)/2
            }
          }
        }
        ###########################################################
        # the last peak is more on the left
      }else if(last_peak<0.2){ 
        # we position the gate on the right of the last peak
        hladr.gate <- deGate(rot$data,fluorochrome.chans['HLA-DR'],use.upper = T,upper = T,alpha = 0.9)
        if(hladr.gate>0.5){ # but the gate is too much on the right
          if(last_peak<0){
            hladr.gate <- (hladr.gate+last_peak)/2
            if((hladr.gate-last_peak)<0.1){
              hladr.gate<-hladr.gate+0.15
            }
          }else{
            hladr.gate <- last_peak*2
          }
        }
      }
    }else{
      ##################### only one peak #######################
      if(all.peaks$Peaks>0.3){ # the peak is on the extreme right,the gate must be on the left
        hladr.gate <- deGate(rot$data,fluorochrome.chans['HLA-DR'], upper=F,alpha = 0.8)
        print(paste0("hladr.gate:",hladr.gate))
        if(hladr.gate<(-0.1)){ # the gate is too much on the left, it must still be on left of the peak but not too much
          hladr.gate <- (all.peaks$Peaks+hladr.gate)/2
          print(all.peaks$Peaks-hladr.gate)
          if((all.peaks$Peaks-hladr.gate)>0.4){
            hladr.gate<-all.peaks$Peaks/2
            if(last_HLADR_peak>2.5){ # the HLDR peak is too much on the right so the hldr gate must be more on the right in these cases
              hladr.gate<-(hladr.gate+all.peaks$Peaks)/2
            }
          }
        }
      }else{ # the peak is on the left,the gate must be on the right
        hladr.gate <- deGate(rot$data,fluorochrome.chans['HLA-DR'], upper=T,alpha = 0.5,tinypeak.removal = 0.3)
        print(paste0("hladr.gate_1:",hladr.gate))
        #show(autoplot(getflowFrame(temp),"eF605-A"))
        if((hladr.gate>0.4)&&(last_HLADR_peak<2.5)){ # gate too much on the right
          hladr.gate_1<-deGate(rot$data,fluorochrome.chans['HLA-DR'], use.percentile = T,percentile = 0.70,tinypeak.removal = 0.3)
          print(hladr.gate_1)
          print(paste0("diff:",hladr.gate_1-all.peaks$Peaks))
          if(last_HLADR_peak>2.4){# the peak of the DRn14n pop is too right respect this peak
            hladr.gate<-(hladr.gate_1+hladr.gate)/2
          }else{
            if((hladr.gate_1-all.peaks$Peaks)>0.6){ # gate too far from the peak
              hladr.gate<-(hladr.gate_1+all.peaks$Peaks)/2
            }else if((hladr.gate_1-all.peaks$Peaks)<=0.3){
              hladr.gate<-deGate(rot$data,fluorochrome.chans['HLA-DR'], use.percentile = T,percentile = 0.80,tinypeak.removal = 0.3)
            }else{
              hladr.gate<-hladr.gate_1
              if(hladr.gate>0 && hladr.gate<0.2){
                hladr.gate<-hladr.gate*1.15
              }
            }
          } 
        }else if(hladr.gate<(-0.4)){ # gate too on the left
          hladr.gate <- deGate(rot$data,fluorochrome.chans['HLA-DR'], upper=T,alpha = 0.1,tinypeak.removal = 0.01) # changed tinypeak thres
          print(paste0("hladr.gate:",hladr.gate))
        } 
        if(last_HLADR_peak>2.6 && hladr.gate<0.1){ # gate must on the right of the peak below (the negaive population below)
          hladr.gate <- deGate(rot$data,fluorochrome.chans['HLA-DR'], upper=T,alpha = 0.001,tinypeak.removal = 0.3)
        }
        
      }
    }
    
    print(paste0("hladr.gate:",hladr.gate))

    monocytes <-  rotate.fd(flowDensity(rot$data, channels = c(fluorochrome.chans['HLA-DR'], fluorochrome.chans['CD14']), position = c(F, F), gates = c(hladr.gate, NA),upper=c(NA,T),use.upper=c(F,T)),angle = rot$theta)
    
    # check NA in monocytes cell population object
    CD14.gatemonoctyes<-deGate(rot$data,fluorochrome.chans['CD14'], use.upper = T,upper = T)
    rot_monocytes_test<-flowDensity(rot$data, channels = c(fluorochrome.chans['HLA-DR'], fluorochrome.chans['CD14']), position = c(F, F), gates = c(hladr.gate,CD14.gatemonoctyes))
    m<-exprs(rot_monocytes_test@flow.frame)

    # ------------- get pop HLADR+CD14-

    #show(autoplot(rot$data,"eF605-A","V500-A"))
    DRp.14n <-  rotate.fd(flowDensity(rot$data, channels = c(fluorochrome.chans['HLA-DR'], fluorochrome.chans['CD14']), position = c(T, NA), gates = c(hladr.gate, NA)),angle = rot$theta)

    poly1 <- SpatialPolygons(list(Polygons(list(Polygon(DRp.14n@filter)), ID=c("c"))))                       
    poly2 <- SpatialPolygons(list(Polygons(list(Polygon(DRn.14n@filter)), ID=c("c"))))
    if(gIntersects(poly1,poly2)){
      n<-colnames(DRp.14n@filter) # nota che i nomi delle colonne di DRp.14n e DRn.14n sono le stesse
      diff.1<-rgeos::gDifference(poly2,poly1)
      DRn.14n@filter<- diff.1@polygons[[1]]@Polygons[[1]]@coords
      colnames( DRn.14n@filter)<-n
    }
    poly1 <- SpatialPolygons(list(Polygons(list(Polygon(monocytes@filter)), ID=c("c"))))
    poly2 <- SpatialPolygons(list(Polygons(list(Polygon(DRn.14n@filter)), ID=c("c"))))
    if (gIntersects(poly1,poly2))
    {
      n<-colnames(DRn.14n@filter)
      diff.1<-rgeos::gDifference(poly1,poly2)
      monocytes@filter<- diff.1@polygons[[1]]@Polygons[[1]]@coords
      colnames(monocytes@filter)<-n
    }

    monocytes<-  flowDensity(fs.66n45p[[i]], channels = c(fluorochrome.chans['HLA-DR'], fluorochrome.chans['CD14']), position = c(T, F),filter=monocytes@filter)
    DRp.14n <-  flowDensity(fs.66n45p[[i]], channels = c(fluorochrome.chans['HLA-DR'], fluorochrome.chans['CD14']), position = c(T, F),filter=DRp.14n@filter)
    DRn.14n <-  flowDensity(fs.66n45p[[i]], channels = c(fluorochrome.chans['HLA-DR'], fluorochrome.chans['CD14']), position = c(T, F),filter=DRn.14n@filter)
    
    # ----- generate polygon gate objects
    mono.poly <- polygonGate(filterId = "HLADR+ CD14+ Monocytes",.gate=monocytes@filter)
    DRp.14n.poly <- polygonGate(filterId = "HLADR+ CD14-",.gate  =DRp.14n@filter)
    DRn.14n.poly <- polygonGate(filterId = "HLADR- CD14-",.gate  =DRn.14n@filter)
    
    return(list(mono.poly=mono.poly,DRp.14n.poly=DRp.14n.poly,DRn.14n.poly=DRn.14n.poly,monocytes=monocytes,DRp.14n=DRp.14n,DRn.14n=DRn.14n))
  })
  list_mono.poly<-lapply(1:length(gs),function(i){
    lista_poly_mono_drp[[i]]$mono.poly
  })
  list_DRp.14n.poly<-lapply(1:length(gs),function(i){
    lista_poly_mono_drp[[i]]$DRp.14n.poly
  })
  list_DRn.14n.poly<-lapply(1:length(gs),function(i){
    lista_poly_mono_drp[[i]]$DRn.14n.poly
  })
  list_monocytes<-lapply(1:length(gs),function(i){
    lista_poly_mono_drp[[i]]$monocytes
  })
  list_DRp_14n<-lapply(1:length(gs),function(i){
    lista_poly_mono_drp[[i]]$DRp.14n
  })
  list_DRn_14n<-lapply(1:length(gs),function(i){
    lista_poly_mono_drp[[i]]$DRn.14n
  })
  
  
  
  names(list_mono.poly)<-sampleNames(gs)
  nodeID0<-add(gs,list_mono.poly,parent="CD45+CD66- (Non Granulocytes)")
  names(list_DRp.14n.poly)<-sampleNames(gs)
  nodeID0<-add(gs,list_DRp.14n.poly,parent="CD45+CD66- (Non Granulocytes)")
  names(list_DRn.14n.poly)<-sampleNames(gs)
  nodeID0<-add(gs,list_DRn.14n.poly,parent="CD45+CD66- (Non Granulocytes)")
  names(list_monocytes)<-sampleNames(gs)
  recompute(gs)
  fs.mono<-getData(gs,"HLADR+ CD14+ Monocytes")
  fs.DRp.14n<-getData(gs,"HLADR+ CD14-")
  fs.DRn.14n<-getData(gs,"HLADR- CD14-")
  return(list(fs.mono=fs.mono,fs.DRp.14n=fs.DRp.14n,fs.DRn.14n=fs.DRn.14n,gs=gs,list_monocytes=list_monocytes,list_DRp_14n=list_DRp_14n,list_DRn_14n=list_DRn_14n))
}


#------  Gating Classical monocytes(CD16-) from monocytes(HLADR+ CD14+) --------------
gating_mono_to_classic_mono<-function(fs.mono,gs,fluorochrome.chans){
  #show(autoplot(fs.66n45p[[1]],"V500-A","FITC-A"))
  set.seed(40)
  list_classicalmono_poly<-lapply(1:length(gs),function(i){
    monocytes<-fs.mono[[i]]
    #CD16monocyte.gate <- deGate(monocytes.flowD, channel = c(fluorochrome.chans['CD16']), use.upper = T, upper = T, tinypeak.removal = 0.2)
    CD16monocyte.gate <- deGate(monocytes, channel = c(fluorochrome.chans['CD16']),upper = T,alpha = 0.3,tinypeak.removal = 0.1)
    print(paste0("CD16monocyte.gate:",CD16monocyte.gate))
    all.peaks<-getPeaks(monocytes, channel = c(fluorochrome.chans['CD16']),tinypeak.removal = 0.1)
    print(all.peaks)
    if(CD16monocyte.gate<=1.7){ # gate too low
      CD16monocyte.gate <- deGate(monocytes, channel = c(fluorochrome.chans['CD16']),use.upper=T,upper = T,tinypeak.removal = 0.1,alpha = 0.05)
      print(CD16monocyte.gate)
      if(CD16monocyte.gate<=1.7){ # gate still too low
        CD16monocyte.gate <- deGate(monocytes, channel = c(fluorochrome.chans['CD16']),use.upper=T,upper = T,tinypeak.removal = 0.1,alpha = 0.001)
        print(CD16monocyte.gate)
      }
      if(CD16monocyte.gate>3.0){ # gate too high
        CD16monocyte.gate <- deGate(monocytes, channel = c(fluorochrome.chans['CD16']),use.upper=T,upper = T,tinypeak.removal = 0.1,alpha = 0.3)
        print(CD16monocyte.gate)
        
      }
      if(CD16monocyte.gate>2.6){ # gate too high
        CD16monocyte.gate <- deGate(monocytes, channel = c(fluorochrome.chans['CD16']),use.upper=T,upper = T,tinypeak.removal = 0.2,alpha = 0.1)
        print(CD16monocyte.gate)
        
      }
      
    }
    if(CD16monocyte.gate>2.4){ # gate too high
      if(max(all.peaks$Peaks)>2.0){ # max peak in high position
        # gate lower
        CD16monocyte.gate <- deGate(monocytes, channel = c(fluorochrome.chans['CD16']),upper = F,tinypeak.removal = 0.1,alpha = 0.3)
      }else if(max(all.peaks$Peaks)<2.0){ # max peak in low position
        # gate higher
        CD16monocyte.gate <- deGate(monocytes, channel = c(fluorochrome.chans['CD16']),use.upper=T,upper = T,tinypeak.removal = 0.1,alpha = 0.01)
      }
      print(CD16monocyte.gate)
    }
    if(CD16monocyte.gate>2.3){
      CD16monocyte.gate<-CD16monocyte.gate*0.80
    }else if(CD16monocyte.gate<1.7){
      CD16monocyte.gate<-CD16monocyte.gate*1.20
    }
    print(paste0("CD16monocyte.gate:",CD16monocyte.gate))
    
    classicalmono<- flowDensity(monocytes, channels = c(fluorochrome.chans['CD14'], fluorochrome.chans['CD16']), position = c(NA, F), gates = c(NA, CD16monocyte.gate))
    nonclassicalmono<- flowDensity(monocytes, channels = c(fluorochrome.chans['CD14'], fluorochrome.chans['CD16']), position = c(NA, T), gates = c(NA, CD16monocyte.gate))
    
    classicalmono.poly <- polygonGate(filterId = "Classical Monocytes",.gate  =classicalmono@filter)
    nonclassicalmono.poly <- polygonGate(filterId = "Non classical Monocytes",.gate  =nonclassicalmono@filter)
    return(list(classicalmono.poly=classicalmono.poly,
                CD16monocyte.gate=CD16monocyte.gate,
                nonclassicalmono.poly=nonclassicalmono.poly,classicalmono=classicalmono,nonclassicalmono=nonclassicalmono))
  })
  list_classicalmono_only_poly <- lapply(1:length(gs),function(i){
    list_classicalmono_poly[[i]]$classicalmono.poly
  })
  list_CD16mon_gate <- lapply(1:length(gs),function(i){
    list_classicalmono_poly[[i]]$CD16monocyte.gate
  })
  list_nonclassicalmono_only_poly <- lapply(1:length(gs),function(i){
    list_classicalmono_poly[[i]]$nonclassicalmono.poly
  })
  list_classicalmono <- lapply(1:length(gs),function(i){
    list_classicalmono_poly[[i]]$classicalmono
  })
  list_nonclassicalmono <- lapply(1:length(gs),function(i){
    list_classicalmono_poly[[i]]$nonclassicalmono
  })
  names(list_classicalmono_only_poly)<-sampleNames(gs)
  nodeID0<-add(gs,list_classicalmono_only_poly,parent="HLADR+ CD14+ Monocytes")
  names(list_nonclassicalmono_only_poly)<-sampleNames(gs)
  nodeID0<-add(gs,list_nonclassicalmono_only_poly,parent="HLADR+ CD14+ Monocytes")
  recompute(gs)
  fs.classmono<- getData(gs,"Classical Monocytes")
  return(list(fs.classmono=fs.classmono,gs=gs,list_CD16mon_gate=list_CD16mon_gate,list_classicalmono=list_classicalmono,list_nonclassicalmono=list_nonclassicalmono))
}  


#---------------------------- Gating HLADR+CD14- to Bcells,pDC,mDC -----------------------------------------------------
gating_drp_14n_to_B_pdc_mdc<-function(fs.drp.14n,gs,fluorochrome.chans,n_cores){
  print("Gating B_pdc_mdc")
  set.seed(40)
  list_poly_B_pdc_mdc<-lapply(1:length(gs),function(i){
    print(identifier(gs[[1]]@data[[i]]))
    #show(autoplot(fs.drp.14n[[1]],"APC-A"))
    #CD11c.gate <- deGate(fs.drp.14n[[i]], channel = c(fluorochrome.chans['CD11c']),bimodal=T) 
    #print(paste0("CD11c.gate:",CD11c.gate))
    all.peaks<-getPeaks(fs.drp.14n[[i]], channel = c(fluorochrome.chans['CD11c']),tinypeak.removal = 0.01)
    print(all.peaks)
    inds<-which(all.peaks$Peaks>0)
    all.peaks$Peaks<-all.peaks$Peaks[inds]
    n_peaks<-length(all.peaks$Peaks)
    CD11c.gate<-sum(all.peaks$Peaks)/n_peaks*1.05
    if(CD11c.gate<2){
      CD11c.gate<-CD11c.gate*1.15
    }
    print(paste0("CD11c.gate:",CD11c.gate))
    # if still incorrect
    if(CD11c.gate<2){
      CD11c.gate<-deGate(fs.drp.14n[[i]], channel = c(fluorochrome.chans['CD11c']),upper = T) 
    }
    # if still incorrect
    if(CD11c.gate<2){
      CD11c.gate<-2.1 # average correct gate across all samples
    }
    print(paste0("CD11c.gate:",CD11c.gate))
    temp <- flowDensity(fs.drp.14n[[i]], channels = c(fluorochrome.chans['CD123'], fluorochrome.chans['CD11c']), position = c(NA, F), gates = c(NA, CD11c.gate))
    temp_2<- flowDensity(fs.drp.14n[[i]], channels = c(fluorochrome.chans['CD123'], fluorochrome.chans['CD11c']), position = c(NA, T), gates = c(NA, CD11c.gate))
    #show(autoplot(getflowFrame(temp_2),"Pe-Cy7-A"))
    #show(autoplot(getflowFrame(temp),"Pe-Cy7-A","APC-A"))
    
    CD123.high.gate <- deGate(temp_2, channel = c(fluorochrome.chans['CD123']), use.upper = T, upper = T, alpha = 0.00001)
    print(paste0("CD123.high.gate:",CD123.high.gate))
    
    CD123.lo.gate <- deGate(temp, channel = c(fluorochrome.chans['CD123']), upper = T,alpha=0.1,tinypeak.removal = 0.1)*1.07
    print(paste0("CD123.lo.gate:",CD123.lo.gate))
    if(CD123.lo.gate>=3.3){
      CD123.lo.gate <- deGate(temp, channel = c(fluorochrome.chans['CD123']),upper = T,alpha = 0.8,tinypeak.removal = 0.1)
    }
    print(paste0("CD123.lo.gate:",CD123.lo.gate))
    # if the gate is still not correct
    if((CD123.lo.gate<2.6) || (CD123.lo.gate>=3.3)){
      all.peaks<-getPeaks(temp, channel = c(fluorochrome.chans['CD123']),tinypeak.removal = 0.1)
      print(all.peaks)
      if((min(all.peaks$Peaks)<2.0) && max(all.peaks$Peaks)>2.5 && length(all.peaks$Peaks)>2){
        inds<-which(all.peaks$Peaks>=1.0)
        all.peaks$Peaks<-all.peaks$Peaks[inds]
        CD123.lo.gate <- (min(all.peaks$Peaks)+max(all.peaks$Peaks))/2
        print(CD123.lo.gate)
      }else if(min(all.peaks$Peaks)>2.0){
        CD123.lo.gate<-min(all.peaks$Peaks)*1.30
        print(CD123.lo.gate)
      }else if((length(all.peaks$Peaks)>2)&&(all.peaks$Peaks[2]>1.8)){
        CD123.lo.gate<-(CD123.lo.gate+max(all.peaks$Peaks))/2
        print(CD123.lo.gate)
      }
      if(CD123.lo.gate<2.5){
        CD123.lo.gate <- deGate(temp, channel = c(fluorochrome.chans['CD123']),use.upper=T,upper = T,alpha = 0.1,tinypeak.removal = 0.1)
        if(CD123.lo.gate<2.1){
          CD123.lo.gate<-CD123.lo.gate*1.40
        }
        if(CD123.lo.gate<2.5){
          CD123.lo.gate<-CD123.lo.gate*1.25
        }
      }
      if(CD123.lo.gate>3.3){
        all.peaks<-getPeaks(temp, channel = c(fluorochrome.chans['CD123']),tinypeak.removal = 0.1)
        if(max(all.peaks$Peaks)>3){
          CD123.lo.gate<-max(all.peaks$Peaks)*0.90
        }else{
          CD123.lo.gate<-max(all.peaks$Peaks)*1.20
        }
      }
    }
    print(paste0("CD123.lo.gate:",CD123.lo.gate))
    mDC <- flowDensity(temp_2, channels = c(fluorochrome.chans['CD123'], fluorochrome.chans['CD11c']), position = c(F, NA), gates = c(CD123.high.gate, NA))
    
    pDC<- flowDensity(temp, channels = c(fluorochrome.chans['CD123'], fluorochrome.chans['CD11c']), position = c(T, F), gates = c(CD123.lo.gate, CD11c.gate))
    
    
    Bcells<- flowDensity(temp, channels = c(fluorochrome.chans['CD123'], fluorochrome.chans['CD11c']), position = c(F, F), gates = c(CD123.lo.gate, CD11c.gate))
    mdc.poly <- polygonGate(filterId = "mDC",.gate  =mDC@filter)
    pdc.poly <- polygonGate(filterId = "pDC",.gate  =pDC@filter)
    bcell.poly <- polygonGate(filterId = "B cells",.gate  =Bcells@filter)
    return(list(mdc.poly=mdc.poly,pdc.poly=pdc.poly,bcell.poly=bcell.poly,
                Bcells=Bcells,
                pDC=pDC,mDC=mDC))
  })
  lista_mdc_poly <- lapply(1:length(list_poly_B_pdc_mdc),function(i){
    list_poly_B_pdc_mdc[[i]]$mdc.poly
  })
  lista_pdc_poly <- lapply(1:length(list_poly_B_pdc_mdc),function(i){
    list_poly_B_pdc_mdc[[i]]$pdc.poly
  })
  lista_bcell_poly <- lapply(1:length(list_poly_B_pdc_mdc),function(i){
    list_poly_B_pdc_mdc[[i]]$bcell.poly
  })
  list_Bcells<-lapply(1:length(list_poly_B_pdc_mdc),function(i){
    list_poly_B_pdc_mdc[[i]]$Bcells
  })
  list_pDC<-lapply(1:length(list_poly_B_pdc_mdc),function(i){
    list_poly_B_pdc_mdc[[i]]$pDC
  })
  list_mDC<-lapply(1:length(list_poly_B_pdc_mdc),function(i){
    list_poly_B_pdc_mdc[[i]]$mDC
  })
  
  
  names(lista_mdc_poly)<-sampleNames(gs)
  nodeID0<-add(gs,lista_mdc_poly,parent="HLADR+ CD14-")
  names(lista_pdc_poly)<-sampleNames(gs)
  nodeID0<-add(gs,lista_pdc_poly,parent="HLADR+ CD14-")
  names(lista_bcell_poly)<-sampleNames(gs)
  nodeID0<-add(gs,lista_bcell_poly,parent="HLADR+ CD14-")
  recompute(gs)
  return(list(gs=gs,lista_mdc_poly=lista_mdc_poly,
              lista_pdc_poly=lista_pdc_poly,lista_bcell_poly=lista_bcell_poly,
              list_Bcells=list_Bcells,
              list_mDC=list_mDC,
              list_pDC=list_pDC))
}



#-------------- Gating HLADR-CD14- to get T cells(CD3+), gdT(CD3+gd+) e CD3-  --------------------------------------------------------
gating_drn_14n_to_T_gdt_cd3<-function(fs.drn.14n,gs,fluorochrome.chans,scat.chans){
  print("gating_tcells_gdt_cd3neg")
  set.seed(40)
  list_t_gdt_cd3_poly<-lapply(1:length(gs), function(i){
    print(identifier(gs[[1]]@data[[i]]))
    CD3.gate <- deGate(fs.drn.14n[[i]], channel = c(fluorochrome.chans['CD3']),upper=F,tinypeak.removal = 0.1,alpha = 0.5)
    print(paste0("CD3.gate:",CD3.gate))
    all.peaks<-getPeaks(fs.drn.14n[[i]], channel = c(fluorochrome.chans['CD3']), tinypeak.removal = 0.1)
    if(CD3.gate>1.8 && CD3.gate<2.5){ # too in central position
      CD3.gate<-(CD3.gate+max(all.peaks$Peaks))/2
      if(CD3.gate>2.8){
        CD3.gate<-CD3.gate*0.95
      }
    }else if(CD3.gate<1.8 || CD3.gate>2.8){ # too much on the left or on the right
      CD3.gate <- deGate(fs.drn.14n[[i]], channel = c(fluorochrome.chans['CD3']),upper=F,tinypeak.removal = 0.01,alpha=0.5)
      print(CD3.gate)
      if(CD3.gate<1.8){ # gate too low
        CD3.gate<-(max(all.peaks$Peaks)+CD3.gate)/2
        print(CD3.gate)
      }
    }
    print(paste0("CD3.gate:",CD3.gate))
    CD3Tcells<- flowDensity(fs.drn.14n[[i]], channels = c(fluorochrome.chans['CD3'], scat.chans['SSC-A']), position = c(T, NA), gates = c(CD3.gate, NA))
    CD3neg <- flowDensity(fs.drn.14n[[i]], channels = c(fluorochrome.chans['CD3'], scat.chans['SSC-A']), position = c(F, NA), gates = c(CD3.gate, NA))
    gd.gate <-deGate(CD3Tcells@flow.frame,fluorochrome.chans['gd'],use.upper=T,upper=T,tinypeak.removal=0.99,alpha=0.05)
    print(paste0("gd.gate:",gd.gate))
    if(gd.gate<2.3){
      gd.gate <- gd.gate*1.10
    }
    if(gd.gate>3.10){
      gd.gate <-deGate(CD3Tcells@flow.frame,fluorochrome.chans['gd'],use.upper=T,upper=T,tinypeak.removal=0.99,alpha=0.11)
      if(gd.gate>3.10){
        gd.gate<-gd.gate*0.90
      }
    }
    print(paste0("gd.gate:",gd.gate))
    temp1<- flowDensity(CD3Tcells@flow.frame, c(fluorochrome.chans['CD3'], fluorochrome.chans['gd']),position = c(NA,T),gates=c(NA,gd.gate))
    temp2<-flowDensity(CD3Tcells@flow.frame, c(fluorochrome.chans['CD3'], fluorochrome.chans['gd']),position = c(NA,F),gates=c(NA,gd.gate))
    cd3.gate.lo.temp1 <- deGate(getflowFrame(temp1),fluorochrome.chans['CD3'],upper=F,use.upper=T,alpha=0.0000001)
    cd3.gate.lo.temp2 <- deGate(getflowFrame(temp2),fluorochrome.chans['CD3'],upper=F,use.upper=T,alpha=0.0000001)
    
    print(paste0("cd3.gate.lo.temp1:",cd3.gate.lo.temp1))
    print(paste0("cd3.gate.lo.temp2:",cd3.gate.lo.temp2))
    gdTcells<-flowDensity(CD3Tcells@flow.frame, c(fluorochrome.chans['CD3'], fluorochrome.chans['gd']),position = c(T,T),gates=c(cd3.gate.lo.temp1,gd.gate))
    gdnegcells<-flowDensity(CD3Tcells@flow.frame, c(fluorochrome.chans['CD3'], fluorochrome.chans['gd']),position = c(T,F),gates=c(cd3.gate.lo.temp2,gd.gate*0.98))
    
    #-----  creiamo i polygon gate objects
    tcell.poly <- polygonGate(filterId = "CD3+ T cells",.gate  = CD3Tcells@filter)
    cd3n.poly <- polygonGate(filterId = "CD3-",.gate  =CD3neg@filter)
    gd.poly <- polygonGate(filterId = "gd T cells",.gate  =gdTcells@filter)
    gdneg.poly<-polygonGate(filterId = "gd- T cells",.gate  =gdnegcells@filter)
    return(list(tcell.poly=tcell.poly,cd3n.poly=cd3n.poly,gd.poly=gd.poly,gdneg.poly=gdneg.poly,CD3.gate=CD3.gate,
                CD3Tcells=CD3Tcells,gdTcells=gdTcells,gdnegcells=gdnegcells,CD3neg=CD3neg))
  })
  list_t_poly<-lapply(1:length(gs),function(i){
    list_t_gdt_cd3_poly[[i]]$tcell.poly
  })
  list_gdt_poly<-lapply(1:length(gs),function(i){
    list_t_gdt_cd3_poly[[i]]$gd.poly
  })
  list_cd3n_poly<-lapply(1:length(gs),function(i){
    list_t_gdt_cd3_poly[[i]]$cd3n.poly
  })
  list_CD3_gate<-lapply(1:length(gs),function(i){
    list_t_gdt_cd3_poly[[i]]$CD3.gate
  })
  list_gdneg_poly<-lapply(1:length(gs),function(i){
    list_t_gdt_cd3_poly[[i]]$gdneg.poly
  })
  list_CD3Tcells<-lapply(1:length(gs),function(i){
    list_t_gdt_cd3_poly[[i]]$CD3Tcells
  })
  list_gdnegcells<-lapply(1:length(gs),function(i){
    list_t_gdt_cd3_poly[[i]]$gdnegcells
  })
  list_gdTcells<-lapply(1:length(gs),function(i){
    list_t_gdt_cd3_poly[[i]]$gdTcells
  })
  list_CD3neg<-lapply(1:length(gs),function(i){
    list_t_gdt_cd3_poly[[i]]$CD3neg
  })
  #----- add gates to the gating tree
  names(list_t_poly)<-sampleNames(gs)
  nodeID0<-add(gs,list_t_poly,parent="HLADR- CD14-")
  recompute(gs)
  names(list_cd3n_poly)<-sampleNames(gs)
  nodeID0<-add(gs,list_cd3n_poly,parent="HLADR- CD14-")
  recompute(gs)
  names(list_gdt_poly)<-sampleNames(gs)
  nodeID0<-add(gs,list_gdt_poly,parent="CD3+ T cells")
  recompute(gs)
  # adding of gd neg population 
  names(list_gdneg_poly)<-sampleNames(gs)
  
  names(list_gdnegcells)<-sampleNames(gs)
  names(list_CD3Tcells)<-sampleNames(gs)
  names(list_gdTcells)<-sampleNames(gs)
  names(list_CD3neg)<-sampleNames(gs)
  
  nodeID0<-add(gs,list_gdneg_poly,parent="CD3+ T cells")
  recompute(gs)
  fs.tcell<- getData(gs,"CD3+ T cells")
  fs.cd3n<-getData(gs,"CD3-")
  return(list(fs.tcell=fs.tcell,fs.cd3n=fs.cd3n,gs=gs,list_gdt_poly=list_gdt_poly,list_CD3_gate=list_CD3_gate,
              list_gdneg_poly=list_gdneg_poly,list_gdnegcells=list_gdnegcells,list_CD3Tcells=list_CD3Tcells,
              list_gdTcells=list_gdTcells,list_CD3neg=list_CD3neg))
  
}

# -------------------------- Gating NKT and NK -------------------------------------------------------------

gating_nkt_nk<-function(fs.tcell,fs.cd3n,gs,fluorochrome.chans){
  print("--- gating nkt")
  set.seed(40)
  #------- get the NKT cells from the CD3+ pop
  list_all_nkt_related_pops<-lapply(1:length(gs),function(i){
    CD16NKT.gate <- deGate(fs.tcell[[i]], channel = c(fluorochrome.chans['CD16']), use.upper = T, upper = T, alpha = 0.05, tinypeak.removal = 0.4)
    if(CD16NKT.gate>3){
      all.peaks<-getPeaks(fs.tcell[[i]], channel = c(fluorochrome.chans['CD16']), tinypeak.removal = 0.01)
      all.cuts <- deGate(fs.tcell[[i]], channel = c(fluorochrome.chans['CD16']), all.cuts = T)
      inds<-which(all.cuts>min(all.peaks$Peaks))
      CD16NKT.gate<-min(all.cuts[inds])
    }
    peaks <-getPeaks(fs.tcell[[i]], channel = c(fluorochrome.chans['CD56']), tinypeak.removal = 0.05)
    #print(peaks$Peaks)
    CD56NKT.gate <-  min(deGate(fs.tcell[[i]], channel = c(fluorochrome.chans['CD56']), tinypeak.removal = 0.05),
                         deGate(fs.tcell[[i]], channel = c(fluorochrome.chans['CD56']), use.upper = T, upper = T, tinypeak.removal = 0.9))
    if(is.na(CD56NKT.gate<peaks$Peaks[which.max(peaks$P.h)])==FALSE){
      if(CD56NKT.gate<peaks$Peaks[which.max(peaks$P.h)]){
        CD56NKT.gate  <- deGate(fs.tcell[[i]], channel = c(fluorochrome.chans['CD56']), use.upper = T, upper = T, tinypeak.removal = 0.9,alpha=.5)
      } 
    }
    #print(CD56NKT.gate)
    if(CD56NKT.gate<1.5){
      CD56NKT.gate  <- deGate(fs.tcell[[i]], channel = c(fluorochrome.chans['CD56']), use.upper = T, upper = T, tinypeak.removal = 0.1,alpha=.5)
      if(CD56NKT.gate>3){
        all.cuts  <- deGate(fs.tcell[[i]], channel = c(fluorochrome.chans['CD56']), all.cuts = T)
        inds<-which(all.cuts>min(peaks$Peaks))
        all.cuts<-all.cuts[inds]
        CD56NKT.gate<-min(all.cuts)
      }
    }
    #print(CD56NKT.gate)
    CD3.16p56p <- flowDensity(fs.tcell[[i]], channels = c(fluorochrome.chans['CD16'], fluorochrome.chans['CD56']), position = c(T, T), gates = c(CD16NKT.gate, CD56NKT.gate))
    CD3.16n56p <- flowDensity(fs.tcell[[i]], channels = c(fluorochrome.chans['CD16'], fluorochrome.chans['CD56']), position = c(F, T), gates = c(CD16NKT.gate, CD56NKT.gate))
    CD3.16n56n <- flowDensity(fs.tcell[[i]], channels = c(fluorochrome.chans['CD16'], fluorochrome.chans['CD56']), position = c(F, F), gates = c(CD16NKT.gate, CD56NKT.gate))
    CD3.16p56n <- flowDensity(fs.tcell[[i]], channels = c(fluorochrome.chans['CD16'], fluorochrome.chans['CD56']), position = c(T, F), gates = c(CD16NKT.gate, CD56NKT.gate))
    
    CD3.16p56p.poly <- polygonGate(filterId = "CD56+CD16+ NKT cells",.gate  =CD3.16p56p@filter)
    CD3.16n56p.poly <- polygonGate(filterId = "CD56-CD16+ NKT cells",.gate  =CD3.16n56p@filter)
    CD3.16p56n.poly <- polygonGate(filterId = "CD56+CD16- NKT cells",.gate  =CD3.16p56n@filter)
    CD3.16n56n.poly <- polygonGate(filterId = "CD56-CD16- NKT cells",.gate  =CD3.16n56n@filter)
    
    return(list(CD3.16p56p=CD3.16p56p,CD3.16n56p=CD3.16n56p,
                CD3.16n56n=CD3.16n56n,CD56NKT.gate=CD56NKT.gate,
                CD16NKT.gate=CD16NKT.gate,
                CD3.16p56p.poly=CD3.16p56p.poly,CD3.16p56n=CD3.16p56n,
                CD3.16n56p.poly=CD3.16n56p.poly,CD3.16p56n.poly=CD3.16p56n.poly,
                CD3.16n56n.poly=CD3.16n56n.poly))
  })
  print("--- gating nk")
  #------- get the NK cells from the CD3- pop
  list_all_nk_pop<-lapply(1:length(gs),function(i){

    f.3n<-fs.cd3n[[i]]
    x<-removeMargins(f.3n,fluorochrome.chans[c('CD16','CD56')],neg=.1);
    CD56NKT.gate <- list_all_nkt_related_pops[[i]]$CD56NKT.gate
    print(paste0("CD56NKT.gate:",CD56NKT.gate))
    temp <- flowDensity(x,fluorochrome.chans[c('CD16','CD56')],position = c(NA,F),gates=c(NA,CD56NKT.gate))
 
    #show(autoplot(getflowFrame(temp),"FITC-A","BV650-A"))

    cd16.gate.hi<-deGate(getflowFrame(temp),fluorochrome.chans['CD16'],upper=T,alpha=0.01,tinypeak.removal = 0.1)
    print(paste0("cd16.gate.hi:",cd16.gate.hi))
    if(cd16.gate.hi>2.5){ # gate too high position
      cd16.gate.hi<-deGate(getflowFrame(temp),fluorochrome.chans['CD16'],upper=T,alpha=0.1,tinypeak.removal = 0.1)
      print(cd16.gate.hi)
      if(cd16.gate.hi>2.5){# still too high position
        cd16.gate.hi<-deGate(getflowFrame(temp),fluorochrome.chans['CD16'],upper=T,alpha=0.5,tinypeak.removal = 0.1)
        print(cd16.gate.hi)
      }
      if(cd16.gate.hi>2.2){
        cd16.gate.hi<-cd16.gate.hi*0.90
      }
      if(cd16.gate.hi>2.5){ # still too high
        all.peaks <- getPeaks(getflowFrame(temp),fluorochrome.chans['CD16'],tinypeak.removal = 0.1)
        cd16.gate.hi<-deGate(getflowFrame(temp),fluorochrome.chans['CD16'],upper=F,alpha=0.1,tinypeak.removal = 0.1)
      }
      print(paste0("cd16.gate.hi.>2.5:",cd16.gate.hi))
    }else if(cd16.gate.hi<1.5){ # gate too low position
      all.peaks <- getPeaks(getflowFrame(temp),fluorochrome.chans['CD16'],tinypeak.removal = 0.1)
      print(all.peaks)
      if(max(all.peaks$Peaks)<1.9){ # last peak on the left
        # gate on the right of the last peak
        cd16.gate.hi<-deGate(getflowFrame(temp),fluorochrome.chans['CD16'],use.upper=T,upper=T,alpha=0.05,tinypeak.removal = 0.1)
        if(cd16.gate.hi>2.6){ # gate too high position
          cd16.gate.hi<-deGate(getflowFrame(temp),fluorochrome.chans['CD16'],use.upper=T,upper=T,alpha=0.3,tinypeak.removal = 0.1)
        }
        if(cd16.gate.hi>2.3){ # still too high
          cd16.gate.hi<-deGate(getflowFrame(temp),fluorochrome.chans['CD16'],use.upper=T,upper=T,alpha=0.5,tinypeak.removal = 0.1)
        }
        if(cd16.gate.hi>2.3){ # still too high
          cd16.gate.hi<-cd16.gate.hi*0.95
        }
        print(cd16.gate.hi)
      }else{ # last peak on the right
        # gate on the left of the last peak
        cd16.gate.hi<-max(all.peaks$Peaks)*0.80
      }
      print(paste0("cd16.gate.hi.<1.5:",cd16.gate.hi))
    }
    print(paste0("cd16.gate.hi:",cd16.gate.hi))
    temp <- flowDensity(x,fluorochrome.chans[c('CD16','CD56')],position = c(NA,T),gates=c(NA,CD56NKT.gate))
    #print(CD56NKT.gate)
    #show(autoplot(getflowFrame(temp),"FITC-A"))
    print(identifier(f.3n))
    all.peaks <- getPeaks(getflowFrame(temp),fluorochrome.chans['CD16'],tinypeak.removal = 0.1)
    print(all.peaks)
    if(length(all.peaks$Peaks)==1){
      if(all.peaks$Peaks<1.9){ # peak on left
        # gate on the right
        cd16.gate.lo <- deGate(getflowFrame(temp),fluorochrome.chans['CD16'],upper=T,alpha=0.05,tinypeak.removal = 0.1)
        if(cd16.gate.lo>2.5){
          cd16.gate.lo <- deGate(getflowFrame(temp),fluorochrome.chans['CD16'],upper=T,alpha=0.3,tinypeak.removal = 0.1)
        }
      }else{ # peak on the right
        # gate on the left
        cd16.gate.lo <- deGate(getflowFrame(temp),fluorochrome.chans['CD16'],upper=F,alpha=0.05,tinypeak.removal = 0.1)
        if(cd16.gate.lo<1.5){ # gate too much on the left
          cd16.gate.lo <- max(all.peaks$Peaks)*0.85
        }
      }
    }else if(length(all.peaks$Peaks)>=2){
      if(max(all.peaks$Peaks)<1.9){ # last peak on the left
        # gate on the right
        cd16.gate.lo <- deGate(getflowFrame(temp),fluorochrome.chans['CD16'],use.upper=T,upper=T,alpha=0.1,tinypeak.removal=0.1)
        if(cd16.gate.lo>2.5){
          cd16.gate.lo <- deGate(getflowFrame(temp),fluorochrome.chans['CD16'],use.upper=T,upper=T,alpha=0.3,tinypeak.removal=0.1)
        }
      }else if(max(all.peaks$Peaks)>=1.9){ # last peak on the right
        # gate on the left
        cd16.gate.lo <- deGate(getflowFrame(temp),fluorochrome.chans['CD16'],upper=T,alpha=0.1,tinypeak.removal = 0.1)
      }
    }
    
    print(paste0("cd16.gate.lo:",cd16.gate.lo))
    quad.2<-flowDensity(f.3n,fluorochrome.chans[c('CD16','CD56')],position = c(T,F),gates=c(cd16.gate.hi,CD56NKT.gate*.9))
    cd56.gate.lo<-deGate(getflowFrame(quad.2),fluorochrome.chans['CD56'],use.upper=T,upper=T,tinypeak.removal = 0.9,count.lim = 5)
    cd56 <- deGate(f.3n,fluorochrome.chans['CD56'],upper=T,tinypeak.removal = 0.0001,bimodal=T)
    if(cd56>=3){
      cd56 <- deGate(f.3n,fluorochrome.chans['CD56'],upper=T,tinypeak.removal = 0.1)
    }
    nk.56n16p <-flowDensity(f.3n,fluorochrome.chans[c('CD16','CD56')],position = c(T,F),gates=c(cd16.gate.hi,cd56))
    nk.56n16n<-flowDensity(f.3n,fluorochrome.chans[c('CD16','CD56')],position = c(F,F),gates=c(cd16.gate.hi,cd56))
    
    NK.temp<-flowDensity(f.3n,fluorochrome.chans[c('CD16','CD56')],position = c(F,T),gates=c(cd16.gate.lo,cd56))
    
    nk.56p16p<-flowDensity(f.3n,fluorochrome.chans[c('CD16','CD56')],position = c(T,T),gates=c(cd16.gate.lo,cd56))
    cd56.gate.sd <- deGate(obj = getflowFrame( NK.temp),channel = fluorochrome.chans['CD56'],percentile=NA,sd.threshold=T,n.sd=1,tinypeak.removal = .8)
   
    cd56.gate.hi <- tail(getPeaks(getflowFrame( NK.temp),fluorochrome.chans['CD56'],tinypeak.removal=.1)$Peaks,1)-.8*sd(exprs(getflowFrame(NK.temp))[,fluorochrome.chans['CD56']])
    print(paste0("cd56.gate.hi:",cd56.gate.hi))
    tryCatch({
      if(cd56.gate.hi<2.8){
        cd56.gate.hi<-tail(getPeaks(getflowFrame( NK.temp),fluorochrome.chans['CD56'],tinypeak.removal=.1)$Peaks,1)*1.10
      }
    },error = function(e) {
      cd56.gate.hi <- tail(getPeaks(getflowFrame( NK.temp),fluorochrome.chans['CD56'],tinypeak.removal=.1)$Peaks,1)*1.10
    })
    # print(identifier(f.3n))
    
    cd56.gate.hi <- max(cd56.gate.sd,cd56.gate.hi)
    nk.56dim16n<- flowDensity( NK.temp,fluorochrome.chans[c('CD16','CD56')],position = c(NA,F),gates=c(NA,cd56.gate.hi))
    nk.56hi <-flowDensity( NK.temp,fluorochrome.chans[c('CD16','CD56')],position = c(NA,T),gates=c(NA,cd56.gate.hi))
    
    # make the polygon gates objects
    nk.56hi.poly <- polygonGate(filterId = "CD56 Hi NK",.gate  =nk.56hi@filter)
    nk.56n16n.poly <- polygonGate(filterId = "CD56-CD16- cells",.gate  =nk.56n16n@filter)
    nk.56n16p.poly <- polygonGate(filterId = "CD56-CD16+ NK",.gate  =nk.56n16p@filter)
    nk.56dim16n.poly <- polygonGate(filterId = "CD56 dim CD16- NK",.gate  =nk.56dim16n@filter)
    nk.56p16p.poly <- polygonGate(filterId = "CD56 dim CD16+ NK",.gate  =nk.56p16p@filter)
    
    
    return(list(nk.56hi.poly=nk.56hi.poly,nk.56n16n.poly=nk.56n16n.poly,
                nk.56n16p.poly=nk.56n16p.poly,nk.56dim16n.poly=nk.56dim16n.poly,
                nk.56p16p.poly=nk.56p16p.poly,nk.56n16n=nk.56n16n,
                NK.temp=NK.temp,nk.56hi=nk.56hi,nk.56dim16n=nk.56dim16n,nk.56p16p=nk.56p16p,nk.56n16p=nk.56n16p))
  })
  list_56hi_poly<-lapply(1:length(gs),function(i){
    list_all_nk_pop[[i]]$nk.56hi.poly
  })
  list_56n16n_poly<-lapply(1:length(gs),function(i){
    list_all_nk_pop[[i]]$nk.56n16n.poly
  })
  list_56n16p_poly<-lapply(1:length(gs),function(i){
    list_all_nk_pop[[i]]$nk.56n16p.poly
  })
  list_56dim16n_poly<-lapply(1:length(gs),function(i){
    list_all_nk_pop[[i]]$nk.56dim16n.poly
  })
  
  list_56p16p_poly<-lapply(1:length(gs),function(i){
    list_all_nk_pop[[i]]$nk.56p16p.poly
  })
  list_nk_56n16n<-lapply(1:length(gs),function(i){
    list_all_nk_pop[[i]]$nk.56n16n
  })
  
  list_NK_temp<-lapply(1:length(gs),function(i){
    list_all_nk_pop[[i]]$NK.temp
  })
  
  list_nk_56hi<-lapply(1:length(gs),function(i){
    list_all_nk_pop[[i]]$nk.56hi
  })
  list_cd16p56pnkt_poly<-lapply(1:length(gs),function(i){
    list_all_nkt_related_pops[[i]]$CD3.16p56p.poly
  })
  list_cd16n56pnkt_poly<-lapply(1:length(gs),function(i){
    list_all_nkt_related_pops[[i]]$CD3.16n56p.poly
  })
  list_cd16n56nnkt_poly<-lapply(1:length(gs),function(i){
    list_all_nkt_related_pops[[i]]$CD3.16n56n.poly
  })
  list_cd16p56nnkt_poly<-lapply(1:length(gs),function(i){
    list_all_nkt_related_pops[[i]]$CD3.16p56n.poly
  })
  
  list_nk.56p16p<-lapply(1:length(gs),function(i){
    list_all_nk_pop[[i]]$nk.56p16p
  })
  list_nk.56dim16n<-lapply(1:length(gs),function(i){
    list_all_nk_pop[[i]]$nk.56dim16n
  })
  list_nk.56n16p<-lapply(1:length(gs),function(i){
    list_all_nk_pop[[i]]$nk.56n16p
  })
  list_CD3.16p56p<-lapply(1:length(gs),function(i){
    list_all_nkt_related_pops[[i]]$CD3.16p56p
  })
  list_CD3.16n56p<-lapply(1:length(gs),function(i){
    list_all_nkt_related_pops[[i]]$CD3.16n56p
  })
  list_CD3.16n56n<-lapply(1:length(gs),function(i){
    list_all_nkt_related_pops[[i]]$CD3.16n56n
  })
  list_CD3.16p56n<-lapply(1:length(gs),function(i){
    list_all_nkt_related_pops[[i]]$CD3.16p56n
  })
  #--------- add gated pops NK and NKT to the gating tree
  names(list_56hi_poly)<-sampleNames(gs)
  nodeID0<-add(gs,list_56hi_poly,parent="CD3-")
  
  names(list_56n16n_poly)<-sampleNames(gs)
  nodeID0<-add(gs,list_56n16n_poly,parent="CD3-")
  
  names(list_56n16p_poly)<-sampleNames(gs)
  nodeID0<-add(gs,list_56n16p_poly,parent="CD3-")
  
  names(list_56dim16n_poly)<-sampleNames(gs)
  nodeID0<-add(gs,list_56dim16n_poly,parent="CD3-")
  
  names(list_56p16p_poly)<-sampleNames(gs)
  nodeID0<-add(gs,list_56p16p_poly,parent="CD3-")
  
  names(list_cd16p56pnkt_poly)<-sampleNames(gs)
  nodeID0<-add(gs,list_cd16p56pnkt_poly,parent="CD3+ T cells")
  
  names(list_cd16n56pnkt_poly)<-sampleNames(gs)
  nodeID0<-add(gs,list_cd16n56pnkt_poly,parent="CD3+ T cells")
  
  names(list_cd16n56nnkt_poly)<-sampleNames(gs)
  nodeID0<-add(gs,list_cd16n56nnkt_poly,parent="CD3+ T cells")
  
  names(list_cd16p56nnkt_poly)<-sampleNames(gs)
  nodeID0<-add(gs,list_cd16p56nnkt_poly,parent="CD3+ T cells")
  
  recompute(gs)
  names(list_nk.56p16p)<-sampleNames(gs)
  names(list_nk.56dim16n)<-sampleNames(gs)
  
  fs.nk.56n16n<-getData(gs,"CD56-CD16- cells")
  return(list(fs.nk.56n16n=fs.nk.56n16n,gs=gs,
              list_all_nkt_related_pops=list_all_nkt_related_pops,
              list_nk_56n16n=list_nk_56n16n,list_NK_temp=list_NK_temp,
              list_nk_56hi=list_nk_56hi,list_CD3.16p56p=list_CD3.16p56p,list_CD3.16n56p=list_CD3.16n56p,
              list_CD3.16n56n=list_CD3.16n56n,list_CD3.16p56n=list_CD3.16p56n,list_nk.56n16p=list_nk.56n16p,
              list_nk.56dim16n=list_nk.56dim16n,list_nk.56p16p=list_nk.56p16p))
}



#----------------- Gating Baso---------------------------------------------------
gating_nk_56n16n_to_Baso<- function(fs.nk.56n16n,gs,fluorochrome.chans){
  set.seed(40)
  list_baso_poly<-lapply(1:length(gs),function(i){
    print(identifier(gs[[1]]@data[[i]]))
    temp <- fs.nk.56n16n[[i]]@exprs[, c(fluorochrome.chans['CD123'], fluorochrome.chans['HLA-DR'])]
    #print(nrow(fs.nk.56n16n[[i]]))
    if(nrow(fs.nk.56n16n[[i]])>10)
    {
      # flowPeaks = A fast and automatic clustering to classify the cells into subpopulations based 
      # on finding the peaks from the overall density function generated by K-means.
      
      flowPeaks.Res3 <- flowPeaks(temp)
      #print(flowPeaks.Res3)
      cd123.gate.nk<- deGate(fs.nk.56n16n[[i]],channel =  fluorochrome.chans['CD123'],bimodal=T,upper=T)
      
      #print("33333")
      basophilcluster.id <- which.min((flowPeaks.Res3$peaks$mu[, 1]  - 3)^2 + (flowPeaks.Res3$peaks$mu[, 2]  - 1)^2)
      #print(basophilcluster.id)
      basophils <- fs.nk.56n16n[[i]]
      
      n_points_baso<-length(which(flowPeaks.Res3$peaks.cluster %in% basophilcluster.id))
      print(paste0("n_points_baso:",n_points_baso))
      if(n_points_baso>1){
        basophils@exprs <-  basophils@exprs[which(flowPeaks.Res3$peaks.cluster %in% basophilcluster.id), ]
      }else{
        basophils@exprs <- t(as.matrix(basophils@exprs[which(flowPeaks.Res3$peaks.cluster %in% basophilcluster.id), ]))
      }
      
      print(max(cd123.gate.nk,min(basophils@exprs[,c(fluorochrome.chans['CD123'])])))
      print( max(basophils@exprs[,c(fluorochrome.chans['HLA-DR'])]))
      basophils_1 <- flowDensity(basophils, channels = c(fluorochrome.chans['CD123'], fluorochrome.chans['HLA-DR']), position = c(T, F), 
                                 gates = c(max(cd123.gate.nk,min(basophils@exprs[,c(fluorochrome.chans['CD123'])])), max(basophils@exprs[,c(fluorochrome.chans['HLA-DR'])])), ellip.gate = T, scale = 0.975)
      output<-any(is.infinite(basophils_1@filter))
      if(output==TRUE){
        cd123_gate<-max(cd123.gate.nk,min(basophils@exprs[,c(fluorochrome.chans['CD123'])]))
        hladr_gate<-max(basophils@exprs[,c(fluorochrome.chans['HLA-DR'])])
        all.peaks_cd123<-getPeaks(fs.nk.56n16n[[i]],fluorochrome.chans["CD123"],tinypeak.removal = 0.3)
        max_peak<-max(all.peaks_cd123$Peaks)
        cd123_gate<-(cd123_gate+max_peak)/2
        if(cd123_gate<1.7){
          cd123_gate<-cd123_gate*1.11
        }
        print(paste0("cd123_gate:",cd123_gate))
        all.peaks_hldr<-getPeaks(fs.nk.56n16n[[i]],fluorochrome.chans["HLA-DR"],tinypeak.removal = 0.3)
        max_peak<-max(all.peaks_hldr$Peaks)
        hladr_gate<-(hladr_gate+max_peak)/2
        print(paste0("hldr_gate:",hladr_gate))
        basophils <- flowDensity(basophils, channels = c(fluorochrome.chans['CD123'], fluorochrome.chans['HLA-DR']), position = c(T, F), 
                                 gates = c(cd123_gate,hladr_gate), ellip.gate = T, scale = 0.975)
      }else{
        basophils<-basophils_1
      }
      print(basophils@filter)
      x_values<-basophils@filter[,1]
      inds<-which(x_values<1.5)
      basophils@filter[inds,1]<-1.5
      basophils<- flowDensity(fs.nk.56n16n[[i]], channels = c(fluorochrome.chans['CD123'], fluorochrome.chans['HLA-DR']), position = c(T,T), filter = basophils@filter)
    }else{
      cd123.gate.nk<- deGate(fs.nk.56n16n[[i]], channel = fluorochrome.chans['CD123'],bimodal=T,upper=T)
      basophils <- flowDensity(fs.nk.56n16n[[i]], channels = c(fluorochrome.chans['CD123'], fluorochrome.chans['HLA-DR']), position = c(T,F),percentile=c(NA,.5))
    }
    # construct the basophils polygon object
    baso.poly <- polygonGate(filterId = "Basophils",.gate  =basophils@filter)
    return(list(baso.poly=baso.poly,basophils=basophils))
  })
  # ------- add gates to the gating tree
  list_basophil_poly<-lapply(1:length(gs),function(i){
    list_baso_poly[[i]]$baso.poly
  })
  list_basophils<-lapply(1:length(gs),function(i){
    list_baso_poly[[i]]$basophils
  })
  names(list_basophil_poly)<-sampleNames(gs)
  names(list_basophils)<-sampleNames(gs)
  nodeID0<-add(gs,list_basophil_poly,parent="CD56-CD16- cells")
  recompute(gs)
  return(list(gs=gs,list_basophil_poly=list_basophil_poly,list_basophils=list_basophils))
}



# ---------------------- Gating Granulocytes --------------------------------------
gating_gran_to_mgran_im_gran<-function(fs.gran,gs,fluorochrome.chans){
  print("gating_gran_to_mgran")
  set.seed(40)
  list_all_gran_poly <- lapply(1:length(gs),function(i){
    print(identifier(gs[[1]]@data[[i]]))
    CD16granulo.gate <- deGate(fs.gran[[i]], channel = c(fluorochrome.chans['CD16']), use.upper = T, upper = F,tinypeak.removal = 0.5,alpha = 0.4)
    print(CD16granulo.gate)
    temp.flowD <- flowDensity(fs.gran[[i]], channels = c(fluorochrome.chans['CD11b'], fluorochrome.chans['CD16']), position = c(NA, T), gates = c(NA, CD16granulo.gate))
    CD11bgranulo.gate <- deGate(temp.flowD, channel = c(fluorochrome.chans['CD11b']), use.upper = T, upper = F, alpha = 0.2,tinypeak.removal = 0.5)
    print(paste0("CD11bgranulo.gate:",CD11bgranulo.gate))
    maturegranulo <- flowDensity(fs.gran[[i]], channels = c(fluorochrome.chans['CD11b'], fluorochrome.chans['CD16']), position = c(T, T), gates = c(CD11bgranulo.gate, CD16granulo.gate))
    immaturegranulo <- flowDensity(fs.gran[[i]], channels = c(fluorochrome.chans['CD11b'], fluorochrome.chans['CD16']), position = c(F, F), gates = c(CD11bgranulo.gate, CD16granulo.gate))
    CD11bposCD16neg.granulo<- flowDensity(fs.gran[[i]], channels = c(fluorochrome.chans['CD11b'], fluorochrome.chans['CD16']), position = c(T, F), gates = c(CD11bgranulo.gate, CD16granulo.gate))
    CD11bnegCD16pos.neutro<- flowDensity(fs.gran[[i]], channels = c(fluorochrome.chans['CD11b'], fluorochrome.chans['CD16']), position = c(F, T), gates = c(CD11bgranulo.gate, CD16granulo.gate))
    # ---- construct polygon gate object
    maturegranulo.poly <- polygonGate(filterId = "CD11b+CD16+ Mature Neutrophils",.gate =maturegranulo@filter)
    immaturegranulo.poly <- polygonGate(filterId = "CD11b-CD16-  Immature Neutrophils 1",.gate =immaturegranulo@filter)
    CD11bposCD16neg.granulo.poly <- polygonGate(filterId = "CD11b+CD16- Granulocytes",.gate =CD11bposCD16neg.granulo@filter)
    CD11bnegCD16pos.neutro.poly <- polygonGate(filterId = "CD11b-CD16+  Immature Neutrophils 2",.gate =CD11bnegCD16pos.neutro@filter)
    return(list(maturegranulo.poly=maturegranulo.poly,immaturegranulo.poly=immaturegranulo.poly,
                CD11bposCD16neg.granulo.poly=CD11bposCD16neg.granulo.poly,
                CD11bnegCD16pos.neutro.poly=CD11bnegCD16pos.neutro.poly,
                CD11bgranulo.gate=CD11bgranulo.gate,
                CD16granulo.gate=CD16granulo.gate,maturegranulo=maturegranulo,immaturegranulo=immaturegranulo,
                CD11bposCD16neg.granulo=CD11bposCD16neg.granulo,CD11bnegCD16pos.neutro=CD11bnegCD16pos.neutro))
  })
  # ------ add gates to the gating tree
  list_mature_gran_poly <- lapply(1:length(gs),function(i){
    list_all_gran_poly[[i]]$maturegranulo.poly
  })
  list_immature_gran_poly <- lapply(1:length(gs),function(i){
    list_all_gran_poly[[i]]$immaturegranulo.poly
  })
  list_CD11bposCD16neg_gran_poly <- lapply(1:length(gs),function(i){
    list_all_gran_poly[[i]]$CD11bposCD16neg.granulo.poly
  })
  list_CD11bnegCD16pos_gran_poly <- lapply(1:length(gs),function(i){
    list_all_gran_poly[[i]]$CD11bnegCD16pos.neutro.poly
  })
  
  list_CD16granulo_gate <- lapply(1:length(gs),function(i){
    list_all_gran_poly[[i]]$CD16granulo.gate
  })
  list_CD11bgranulo_gate <- lapply(1:length(gs),function(i){
    list_all_gran_poly[[i]]$CD11bgranulo.gate
  })
  list_maturegranulo <- lapply(1:length(gs),function(i){
    list_all_gran_poly[[i]]$maturegranulo
  })
  list_immaturegranulo <- lapply(1:length(gs),function(i){
    list_all_gran_poly[[i]]$immaturegranulo
  })
  list_CD11bposCD16neg.granulo <- lapply(1:length(gs),function(i){
    list_all_gran_poly[[i]]$CD11bposCD16neg.granulo
  })
  list_CD11bnegCD16pos.neutro <- lapply(1:length(gs),function(i){
    list_all_gran_poly[[i]]$CD11bnegCD16pos.neutro
  })
  names(list_mature_gran_poly)<-sampleNames(gs)
  nodeID0<-add(gs,list_mature_gran_poly,parent="Granulocytes")
  names(list_immature_gran_poly)<-sampleNames(gs)
  nodeID0<-add(gs,list_immature_gran_poly,parent="Granulocytes")
  names(list_CD11bposCD16neg_gran_poly)<-sampleNames(gs)
  nodeID0<-add(gs,list_CD11bposCD16neg_gran_poly,parent="Granulocytes")
  names(list_CD11bnegCD16pos_gran_poly)<-sampleNames(gs)
  nodeID0<-add(gs,list_CD11bnegCD16pos_gran_poly,parent="Granulocytes")
  names(list_maturegranulo)<-sampleNames(gs)
  recompute(gs)
  fs.11p16p<-getData(gs,"CD11b+CD16+ Mature Neutrophils")
  return(list(fs.11p16p=fs.11p16p,gs=gs,list_CD16granulo_gate=list_CD16granulo_gate,
              list_CD11bgranulo_gate=list_CD11bgranulo_gate,
              list_maturegranulo=list_maturegranulo,list_immaturegranulo=list_immaturegranulo,
              list_CD11bposCD16neg.granulo=list_CD11bposCD16neg.granulo,list_CD11bnegCD16pos.neutro=list_CD11bnegCD16pos.neutro))
}



# ---------------------------- Gating CD64+ ----------------------------------------------------------------------
gating_mneutro_to_CD64pos<-function(fs.11p16p,gs,fluorochrome.chans){
  set.seed(40)
  list_all_CD64_poly <-lapply(1:length(gs),function(i){
    print(identifier(gs[[1]]@data[[i]]))
    all.peaks<-getPeaks(fs.11p16p[[i]],fluorochrome.chans["CD64"],tinypeak.removal = 0.3)
    print(all.peaks)
    ind<-which(all.peaks$Peaks==max(all.peaks$Peaks))
    right_peak<-all.peaks$Peaks[ind]
    if(right_peak<1.7){
      cd64.gate<- deGate(fs.11p16p[[i]],fluorochrome.chans["CD64"],upper=T,alpha = 0.8,tinypeak.removal = 0.01)
    }else{
      cd64.gate<- right_peak*0.75
    }
    print(paste0("cd64.gate:",cd64.gate))
    # --- construct polygon gate object
    CD64n.poly <- polygonGate(filterId ="CD64-",.gate = flowDensity(fs.11p16p[[1]],fluorochrome.chans[c("CD64","CD11b")],position = c(F,NA),gates=c(cd64.gate,NA))@filter)
    CD64p.poly <- polygonGate(filterId ="CD64+",.gate = flowDensity(fs.11p16p[[1]],fluorochrome.chans[c("CD64","CD11b")],position = c(T,NA),gates=c(cd64.gate,NA))@filter)
    return(list(CD64n.poly=CD64n.poly,CD64p.poly=CD64p.poly,cd64.gate=cd64.gate))
  })
  list_CD64n_poly<-lapply(1:length(gs),function(i){
    list_all_CD64_poly[[i]]$CD64n.poly
  })
  list_CD64p_poly<-lapply(1:length(gs),function(i){
    list_all_CD64_poly[[i]]$CD64p.poly
  })
  list_CD64_gate<-lapply(1:length(gs),function(i){
    list_all_CD64_poly[[i]]$cd64.gate
  })
  
  #----  add gates to the gating Tree
  names(list_CD64n_poly)<-names(list_CD64p_poly)<-sampleNames(gs)
  nodeID0<-add(gs,list_CD64p_poly,parent="CD11b+CD16+ Mature Neutrophils")
  nodeID0<-add(gs,list_CD64n_poly,parent="CD11b+CD16+ Mature Neutrophils")
  recompute(gs)
  return(list(gs=gs,list_CD64_gate=list_CD64_gate))
}

