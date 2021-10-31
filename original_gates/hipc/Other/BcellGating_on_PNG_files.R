

#########################################################################################################3
set.seed(40)


################################################## Inizio del Gating ######################################################


#----- Gating dei margins----------------
margin_gating <- function(){
  margin.gates <- fsApply(fs,removeMargins,c("FSC-A","SSC-A"),return.gate=TRUE)
  rg <- lapply(1:nrow(margin.gates), function(x) return(rectangleGate(filterId="Margin", "FSC-A"=c(-1, margin.gates[x,1]),
                                                                      "SSC-A"=c(-1, margin.gates[x,2]))))
  # assegnamo come nome di rg(che attualmente e' NULL) il nome dei samples di fs(che c'è ne solo uno quindi un solo nome)
  names(rg) <-sampleNames(gs)
  # Aggiugiamo questo popolazione/nodo/gate al GatingSet(e quindi al gating tree dei suoi GatingHierarchy,in questo caso solo 1)
  node_margin<-add(gs, rg)#,"Margin")
  # add presenta la seguente struttura:
  # add(wf, action, ...) dove:
  # wf = A GatingHierrarchy or GatingSet
  # action = A filter or a list of filters to be added to the GatingHierarchy or GatingSet.
  # Eseguiamo effettivamente il Gating sul flowSet del GatingSet in base alle informazioni che abbiamo appena caricato.
  recompute(gs)
  fs.marg <- getData(gs,"Margin")
  return(list(gs=gs,fs.marg=fs.marg))
}







# #------- Compensation e transformation sul GatingSet  -------------------

comp_and_transform<- function(){
  # ATTENZIONE ricorda che la compensazione nel Gating Set deve essere eseguita con questo codic
  gs<- compensate.flow(gs,comp = comp)
  # NON con compensate(gs,as.matrix(comp)) se no GatingSet2flowJo ti da' errore
  # l'output è il GatingSet con il flowSet compensato
  recompute(gs)
  # Aggiorniamo il GatingSet nella versione con il flowSet compensato
  # usiamo la funzione trasform.flow del flowPrep.R file per eseguire la trasformazione sul GatingSet.
  gs <-transform.flow(gs,remove.outliers=F,trans.chans = NULL)
  # l'output è il GatingSet con il flowSet trasformato.
  recompute(gs)
  fs.marg <- getData(gs,"Margin")
  return(list(gs=gs,fs.marg=fs.marg))
}


#------- Gating non marginal to cleaned ---------------------------
cleaning<-function(){
  Clean <- fsApply(fs.marg, function(f){
    # # Rimuoviamo questo codice percheè l'ho abbiamo eseguito fuori con altre funzioni
    # # eseguiamo la compensazione
    # f1 <- compensate(f,as.matrix(comp))
    # # eseguiamo la trasformazione
    # f1 <- estimate.logicle(fs.raw = as(f1,"flowSet"),med=F,estimate=T,trans.chans = c(7:19))
    # f1 <- f1[[1]]
    f1 <- f
    print(identifier(f1))
    # eseguiamo il cleaning con flowCut,al posto di settare Segment=10000 settiamo Segment=1000,perche'
    # un sample ha dato errore,affermando che il Segment era troppo grande.
    res.tvsf <- flowCut(f1, Segment=10000, Channels=NULL, Directory=paste0(path.output,"/Cleaning/"),
                        FileID= paste0(identifier(f),"_"), Plot="All")
    return(res.tvsf)
  })
  clean.clust <- lapply(1:length(Clean), function(x){
    vec<-rep(0,nrow(fs.marg[[x]]))
    if (length(Clean[[x]]$ind)>0)
    {
      vec[Clean[[x]]$ind]<-1
    }
    vec <- as.factor(vec)
    levels(vec) <- c("0", "1")
    return(vec)
  })
  names(clean.clust)<-sampleNames(gs)
  add(gs,clean.clust,parent = "Margin", name = "Clean")
  recompute(gs)
  fs.clean <- getData(gs,"Clean_0")
  return(list(fs.clean=fs.clean,gs=gs,Clean=Clean))
}





#show(autoplot(fs.clean[[1]],"FSC-A","FSC-H"))

# ---------------------------- Gating from cleaned cells to singlets --------------------------
gating_time_to_singlets<- function(){
  Singlet <- fsApply(fs.clean, function(f) {
    print(identifier(f))
    channels <- c("FSC-A", "FSC-H")
    rot <- rotate.data(f, channels,theta = -pi/4)$data
    #show(autoplot(rot,"FSC-A","FSC-H"))
    #show(autoplot(rot,"FSC-A"))
    gate.2 <- deGate(rot, channels[1],percentile=NA,tinypeak.removal = .98,twin.factor = .2,magnitude = .05, upper=T, alpha=.007,verbose = F)*1.06 
    # mio codice
    #gate.2 <- deGate(rot, channels[1],upper=T,alpha=0.05,twin.factor = .2) # uso twin.factor per rimuovere il picco piccolo nel sample A10.,cosi' il treshold e' piu' a destra di questo picco.
    #print(gate.2)
    singlets<- rotate.fd(flowDensity(rot,channels,position = c(F,NA),gates=c(gate.2,NA),verbose = F),angle = -pi/4)
    return(singlets)
  })
  
  #----- aggiorniamo il gatingSet
  singlets_filter_list <- lapply(Singlet,function(x) x@filter)
  sngl.poly <- lapply(1:length(singlets_filter_list),function(i){
    polygonGate(filterId = "Singlets",.gate=singlets_filter_list[i])
  })
  names(sngl.poly)<-sampleNames(gs)
  node_singlets<-add(gs,sngl.poly,parent="Clean_0")
  recompute(gs)
  fs.sngl <- getData(gs,"Singlets")
  return(list(gs=gs,fs.sngl=fs.sngl,Singlet=Singlet))
}



# ------------------------ Gating from Singlets to Size and Beads----------------------
#  nel caso delle Bcells Gambia,c'e' un problema nel gating dei size,che sfora verso le beads
# visualizziamo dati fs.singlets respect to SSC-A:show(autoplot(fs.sngl[[2]],"SSC-A"))

# entrambi i channels:show(autoplot(fs.sngl[[2]],"FSC-A","SSC-A"))

gating_singlets_to_size_beads <- function(){
  
  size <- fsApply(fs.sngl, function(x){
    print(paste("Gating size of ",identifier(x)))
    
    #x<-removeMargins(x,c("FSC-A","SSC-A"),neg = 500,sens = .999)
    x.rot <- rotate.data(data = x,chans =c("FSC-A","SSC-A"),theta=-pi/2.5)$data
    peaks <-getPeaks(x.rot, "SSC-A",tinypeak.removal = .02)$Peaks
    if(peaks[which.min(peaks-55000)]>65000){
      ss.lo <-  deGate(x.rot, "SSC-A", upper=F,use.upper=T, percentile=NA,alpha=.5)
    }else{
      ss.lo <- deGate(x.rot, "SSC-A", upper=NA, percentile=NA,alpha=0.025,all.cut=T,tinypeak.removal = .02,n.sd=.4)[1]*.95
    }
    temp <-flowDensity(x.rot,channels =c("FSC-A","SSC-A"), position=c(NA,T),gates=c(NA,ss.lo))
    temp@flow.frame <- rotate.data(getflowFrame(temp),chans =c("FSC-A","SSC-A"),theta =pi/2.5)$data
    ss.hi <- deGate(getflowFrame(temp), "SSC-A", upper=T, percentile=NA,alpha=0.1,all.cut=T,tinypeak.removal = 1/40)
    ss.hi <-ss.hi[which.min(abs(ss.hi-170000))]
    size <- flowDensity(temp,channels =c("FSC-A","SSC-A"), position=c(NA,F),gates=c(NA,ss.hi))
    size@proportion<-size@cell.count/nrow(x)*100

    return(size)
  })
  beads <- fsApply(fs.sngl, function(x){
    print(paste("Gating beads of ",identifier(x)))
    temp <-flowDensity(x,channels =c("FSC-A","SSC-A"), position=c(NA,T),gates=c(NA,size[[identifier(x)]]@gates[2]))
    gate <- deGate(temp, "FSC-A", percentile=NA,alpha=.01,use.upper=T,upper=T,tinypeak.removal = .8)
    temp2 <- flowDensity(temp,channels =c("FSC-A","SSC-A"), position=c(F,NA), gates=c(gate,NA))
    bead <- flowDensity(temp2,channels =c("FSC-A","SSC-A"), position=c(F,NA), use.percentile=c(T,F),percentile=c(.999,NA),ellip.gate=T)
    bead@proportion<-bead@cell.count/nrow(x)*100
    return(bead)
  })
  
  # ----- aggiorniamo il GatingSet
  
  size_filter_list <- lapply(size,function(x) x@filter)
  size.poly <- lapply(1:length(size_filter_list),function(i){
    polygonGate(filterId = "Size",.gate=size_filter_list[i])
  })
  names(size.poly)<-sampleNames(gs)
  node_size<-add(gs,size.poly,parent="Singlets")
  recompute(gs)
  # ATTENZIONE: non c'e bisogno di performare di nuovo la compensation ne transformation come ha fatto Mehrnoush perche' io ho compensato  e trasformato prima a partire dal gating Set
  # lei invece ha seguito un altra strada.
  
  beads_filter_list <- lapply(beads,function(x) x@filter)
  beads.poly <- lapply(1:length(beads_filter_list),function(i){
    polygonGate(filterId = "Beads",.gate=beads_filter_list[i])
  })
  names(beads.poly)<-sampleNames(gs)
  node_beads<-add(gs,beads.poly,parent="Singlets")
  recompute(gs)
  fs.size <- getData(gs,"Size")
  return(list(gs=gs,fs.size=fs.size,size=size,beads=beads))
}

# ------------------------ Gating from size to Live cells ----------------------
# usiamo la funzione averageGates per trovare i gate mediani del flowSet. Vuole come input un vettore numerico,infatti 
# fsApply(fs.size,function(x) tail(deGate(x,"APC-eF780-A",upper=T,use.upper = T, tinypeak.removal = .5,alpha=0.05),1)) dà un vettore numerico
# infatti ricorda : head,tail Show first/last elements of the object,in this case (tail,...,1) otteniamo solo l'ultimo elemento,infatti scriviamo 1.
# Dunque fsApply(fs.size,function(x) tail(deGate(x,"APC-eF780-A",upper=T,use.upper = T, tinypeak.removal = .5,alpha=0.05),1)) è una lista di treshold,
# uno per ogni flowFrame
gating_size_to_live<-function(){
  live.gate <- averageGates(fsApply(fs.size,function(x) tail(deGate(x,"APC-eF780-A",upper=T,use.upper = T, tinypeak.removal = .5,alpha=0.05),1)),sd.coeff = 3.5)
  names(live.gate) <- sampleNames(gs)
  # visualizziamo pop live rispetto a viability dye e SSC-A: show(autoplot(fs.size,"APC-eF780-A")),show(autoplot(fs.size[[1]],"SSC-A"))
  # rispetto a entrambi: show(autoplot(fs.size[[1]],"APC-eF780-A","SSC-A"))
  live <- fsApply(fs.size, function(x){
    live <-flowDensity(x,channels =c("APC-eF780-A","SSC-A"), position=c(F,NA),gates=c(live.gate[identifier(x)],NA))
    return(live)
  })
  
  # ----- aggiorniamo il GatingSet
  
  live_filter_list <- lapply(live,function(x) x@filter)
  live.poly <- lapply(1:length(live_filter_list),function(i){
    polygonGate(filterId = "Live",.gate=live_filter_list[i])
  })
  names(live.poly)<-sampleNames(gs)
  node_live<-add(gs,live.poly,parent="Size")
  recompute(gs)
  fs.live <- getData(gs,"Live")
  return(list(gs=gs,fs.live=fs.live,live=live))
}

# ------------------------ Gating from Live to Non Granulocytes----------------------
gating_live_to_non_Granulocytes<- function(){
  # dalla pop madre live dobbiamo ottenere la pop dei non Granulocytes(CD66-CD14-)
  cd66.gate <- averageGates(fsApply(fs.live,function(x) deGate(x,channels.ind[2],upper=T,tinypeak.removal = 0.1)),sd.coeff = 2)
  names(cd66.gate) <- sampleNames(gs)
  # rappresentiamo i fs.live rispetto a CD66 e Cd14: show(autoplot(fs.live[[1]],"BV711-A")), show(autoplot(fs.live[[1]],"V500-A"))
  # rispetto a entrambi: show(autoplot(fs.live[[1]],"BV711-A","V500-A"))
  nonGran <- fsApply(fs.live, function(x){
    
    temp <- exprs(x)[,channels.ind[2:3]]
    fp <- flowPeaks(temp)
    clusterid <- which(fp$peaks$mu[,1]<cd66.gate[identifier(x)]*.95)
    non.g <- x
    non.g@exprs <- non.g@exprs[which(fp$peaks.cluster %in% clusterid), ]
    c.hull <- chull(non.g@exprs[,channels.ind[2:3]])
    c.hull <- non.g@exprs[c(c.hull,c.hull[1]),channels.ind[2:3]]
    colnames(c.hull)<- colnames(x)[channels.ind[2:3]]
    
    cell.population <- new("CellPopulation",
                           flow.frame=non.g,
                           proportion=(nrow(non.g)/nrow(x)*100),
                           cell.count=nrow(non.g),
                           channels=colnames(x)[channels.ind[2:3]],
                           position=c(F,NA),
                           gates=c(cd66.gate[identifier(x)],NA),
                           filter=c.hull,
                           index=which(fp$peaks.cluster %in% clusterid))
    
    return(cell.population)})
  
  # ----- aggiorniamo il GatingSet
  
  nonGran_filter_list <- lapply(nonGran,function(x) x@filter)
  nonGran.poly <- lapply(1:length(nonGran_filter_list),function(i){
    polygonGate(filterId = "NonGran",.gate=nonGran_filter_list[i])
  })
  names(nonGran.poly)<-sampleNames(gs)
  node_nonGran<-add(gs,nonGran.poly,parent="Live")
  recompute(gs)
  fs.ngran <- getData(gs,"NonGran")
  return(list(gs=gs,fs.ngran=fs.ngran,nonGran=nonGran))
  
}




# ------------------------ Gating from non-Granulocytes to Bcells(CD19+)----------------------
# visualizziamo pop non gran rispetto a CD19:show(autoplot(fs.ngran[[2]],"APC-A","SSC-A"))
# show(autoplot(fs.ngran[[3]],"APC-A","BV605-A"))
# show(autoplot(fs.ngran[[1]],"APC-A"))
# rispetto a SSC-A:show(autoplot(fs.ngran[[2]],"SSC-A"))
# proporzioni manual count file A10,A11,A12: 0.05576243,0.04519665,0.06943639
gating_nonGran_to_Bcells<- function(){
  
  ss.ind <- grep("SSC-A",colnames(fs[[1]]))
  cd19.gate <- fsApply(fs.ngran,deGate, channels.ind[6],tinypeak.removal=1/50) # se vuoi prende la pop il piu' a destra possibile imposta alpha=0.1 e togli *0.90
  names(cd19.gate)<-sampleNames(gs)
  Bcell <- fsApply(fs.ngran, function(x){
    cd19<- flowDensity(x,c(channels.ind[6],ss.ind),position = c(T,NA),gates=c(cd19.gate[identifier(x)],NA)) 
    # in pratica al posto di cd19.gate[identifier(x)]*.95 ho tolto il 0.95,per spingere il treshold piu' a destra
    return(cd19)
  })
  # ----- aggiorniamo il GatingSet
  
  Bcell_filter_list <- lapply(Bcell,function(x) x@filter)
  Bcell.poly <- lapply(1:length(Bcell_filter_list),function(i){
    polygonGate(filterId = "Bcells",.gate=Bcell_filter_list[i])
  })
  names(Bcell.poly)<-sampleNames(gs)
  node_Bcell<-add(gs,Bcell.poly,parent="NonGran")
  recompute(gs)
  fs.19 <- getData(gs,"Bcells")
  return(list(gs=gs,fs.19=fs.19,Bcell=Bcell))
}


# ------------------------ Gating from Bcells to PC----------------------
# show(autoplot(fs.19[[3]],"PE-Cy5-A","V450-A"))
# show(autoplot(fs.19[[1]],"PE-Cy5-A"))
# show(autoplot(fs.19[[2]],"V450-A"))
gating_bcells_to_pc<-function(){
  
  PC <- fsApply(fs.19, function(x){
    # counter.line richiama la funzione counterLines di R che calcola le countour lines dato z(densita' bivariata),x e y.
    # di cui la funzione counter.line prevede di prendere le coordinate che si riferiscono alla prima counter line (which.line=1)
    cd138.gate <-max(contour.line(x,channels.ind[7:8])[,2])
    cd38.gate <- deGate(x, channels.ind[7],upper=T, use.upper = T, magnitude = .7,alpha = .5)
    pc <- flowDensity(x,channels.ind[7:8],position = c(T,T),gates=c(cd38.gate,cd138.gate))
    pc@filter <- cbind(c(cd38.gate,max(exprs(x)[,channels.ind[7]],na.rm = T)*.99,max(exprs(x)[,channels.ind[7]],na.rm = T)*.99,cd38.gate,cd38.gate),
                       c(cd138.gate,cd138.gate,max(exprs(x)[,channels.ind[8]],na.rm = T)*.98,max(exprs(x)[,channels.ind[8]],na.rm = T)*.98,cd138.gate))
    colnames(pc@filter)<-colnames(x)[channels.ind[7:8]]
    # print(cd138.gate)
    # print(cd38.gate)
    return(pc)
  })
  
  # ----- aggiorniamo il GatingSet
  
  PC_filter_list <- lapply(PC,function(x) x@filter)
  PC.poly <- lapply(1:length(PC_filter_list),function(i){
    polygonGate(filterId = "PC",.gate=PC_filter_list[i])
  })
  names(PC.poly)<-sampleNames(gs)
  node_PC<-add(gs,PC.poly,parent="Bcells")
  recompute(gs)
  return(list(gs=gs,PC=PC))
}



# ------------------------Gating from Bcells to PB----------------------
# proporzioni manual count A10,A11,A12: 0.004451864,0.006901084,0.003581021
#show(autoplot(fs.19[[3]],"PE-Cy5-A","PE-Cy7-A"))
#show(autoplot(fs.19[[3]],"PE-Cy5-A"))
#show(autoplot(fs.19[[1]],"PE-Cy7-A"))
gating_bcells_to_PB<-function(){
  PB <- fsApply(fs.19, function(x){
    cd38.gate <-deGate(x,channels.ind[7],upper = F, alpha=.9, use.upper=T)
    cd27.gate <- deGate(x, channels.ind[12],upper=T, alpha=.08,tinypeak.removal = 0.9)
    pb <- flowDensity(x,channels.ind[c(7,12)],position = c(T,T),gates=c(cd38.gate,cd27.gate))
    pb@filter <- cbind(c(cd38.gate,max(exprs(x)[,channels.ind[7]],na.rm = T)*.99,max(exprs(x)[,channels.ind[7]],na.rm = T)*.99,cd38.gate,cd38.gate),
                       c(cd27.gate,cd27.gate,max(exprs(x)[,channels.ind[12]],na.rm = T)*.99,max(exprs(x)[,channels.ind[12]],na.rm = T)*.99,cd27.gate))
    colnames(pb@filter)<-colnames(x)[channels.ind[c(7,12)]]
    return(pb)
  })
  
  # ----- aggiorniamo il GatingSet
  PB_filter_list <- lapply(PB,function(x) x@filter)
  PB.poly <- lapply(1:length(PB_filter_list),function(i){
    polygonGate(filterId = "PB",.gate=PB_filter_list[i])
  })
  names(PB.poly)<-sampleNames(gs)
  node_PB<-add(gs,PB.poly,parent="Bcells")
  recompute(gs)
  return(list(gs=gs,PB=PB))
}



# ------------------------Gating from Bcells to CD19/20----------------------
# proportional manual count A10,A11,A12 : 
# rappresentiamo flowset rispetto a cd19 e cd20: show(autoplot(fs.19[[6]],"APC-A")), show(autoplot(fs.19[[3]],"BV786-A"))
# rispetto a entrambi: show(autoplot(fs.19[[1]],"APC-A","BV786-A"))
#show(autoplot(fs.19[[2]],"BV786-A"))

gating_bcells_to_CD20p<-function(){
  CD19.20 <- fsApply(fs.19, function(x){
    cd19.gate <-deGate(x,channels.ind[6],upper = T,use.upper=T, alpha=.05)
    cd20.gate <- deGate(x, channels.ind[9],upper=F, alpha=.01,twin.factor = .8,tinypeak.removal = 1/70)
    cd20.pos<- flowDensity(x,channels.ind[c(6,9)],position = c(F,T),gates=c(cd19.gate,cd20.gate))
    cd20.neg<- flowDensity(x,channels.ind[c(6,9)],position = c(F,F),gates=c(cd19.gate,cd20.gate))
    # pb@filter <- cbind(c(cd38.gate,max(exprs(x)[,channels.ind[7]],na.rm = T)*.99,max(exprs(x)[,channels.ind[7]],na.rm = T)*.99,cd38.gate,cd38.gate),
    #                    c(cd27.gate,cd27.gate,max(exprs(x)[,channels.ind[12]],na.rm = T)*.99,max(exprs(x)[,channels.ind[12]],na.rm = T)*.99,cd27.gate))
    # colnames(pc@filter)<-colnames(x)[channels.ind[c(7,12)]]

    return(list(pos=cd20.pos, neg=cd20.neg))
  })
  
  
  
  
  # ----- aggiorniamo il GatingSet
  
  CD19.20pos_filter_list <- lapply(CD19.20,function(x) x$pos@filter)
  CD19.20pos.poly <- lapply(1:length(CD19.20pos_filter_list),function(i){
    polygonGate(filterId = "CD19.20pos",.gate=CD19.20pos_filter_list[i])
  })
  names(CD19.20pos.poly)<-sampleNames(gs)
  node_CD19.20pos<-add(gs,CD19.20pos.poly,parent="Bcells")
  recompute(gs)
  
  CD19.20neg_filter_list <- lapply(CD19.20,function(x) x$neg@filter)
  CD19.20neg.poly <- lapply(1:length(CD19.20neg_filter_list),function(i){
    polygonGate(filterId = "CD19.20neg",.gate=CD19.20neg_filter_list[i])
  })
  names(CD19.20neg.poly)<-sampleNames(gs)
  node_CD19.20neg<-add(gs,CD19.20neg.poly,parent="Bcells")
  recompute(gs)
  fs.20 <- getData(gs,"CD19.20pos")
  return(list(fs.20=fs.20,gs=gs,CD19.20=CD19.20))
}

# vediamo differenze
# plotDens(fs.19[[1]],c(7,12))
# lines(CD19.20[[1]]$pos@filter,type="l",lwd=2)
# plotDens(getData(gs[[1]],"Bcells"),c(7,12))
# lines(CD19.20[[1]]$pos@filter,type="l",lwd=2)
# lines(CD19.20[[1]]$neg@filter,type="l",lwd=2)
# 
# 
# 
# plotDens(fs.19[[2]],c(7,12))
# lines(CD19.20[[2]]$pos@filter,type="l",lwd=2)
# plotDens(getData(gs[[2]],"Bcells"),c(7,12))
# lines(CD19.20[[2]]$pos@filter,type="l",lwd=2)
# lines(CD19.20[[1]]$neg@filter,type="l",lwd=2)




# ------------------------Gating from Bcells to IgM/IgD----------------------
#show(autoplot(fs.19[[3]],"PE-CF594-A","FITC-A"))
# show(autoplot(fs.19[[2]],"PE-CF594-A"))
# show(autoplot(fs.19[[2]],"FITC-A"))
gating_bcells_to_igm_igd<- function(i){
  # igm.gate <- averageGates(as.numeric(fsApply(fs.19,deGate,channels.ind["IgM"],percentile=NA,upper=F)),sd.coeff = 2)
  # igd.gate <- averageGates(as.numeric(fsApply(fs.19,function(x) return(deGate(x,channels.ind["IgD"],use.upper=T,upper=T,tinypeak.removal = 0.5)[1]))),sd.coeff = 2)
  
  igm.gate <- averageGates(as.numeric(fsApply(fs.19,deGate,channels.ind["IgM"],percentile=NA,upper=F)),sd.coeff = 2)
  igd.gate <- averageGates(as.numeric(fsApply(fs.19,function(x) return(deGate(x,channels.ind["IgD"],after.peak = F,percentile=NA,upper=F)[1]))),sd.coeff = 2)

  
  # igd.gate<-igd.gate*1.15
  # # aggiungiamo del codice per fare in modo che funzioni sia PNG samples che su Gambia samples 
  # all_peaks <- fsApply(fs.19,getPeaks,channels.ind["IgD"])
  # max_indx <- max(which(all_peaks$V1$Peaks>2.1)) # indice del picco più spostato a destra dopo 2.1
  # se l'altezza del picco più a destra è maggiore dell'altezza del primo picco della distribuzione,
  # prendiamo il threshold il più a sinistra possibile e non a destra come nel primo codice
  # if(all_peaks$V1$P.h[max_indx] > all_peaks$V1$P.h[1]){
  #   igd.gate <- averageGates(as.numeric(fsApply(fs.19,function(x) return(deGate(x,channels.ind["IgD"],after.peak = F,upper=F)[1]))),sd.coeff = 2)
  # }

  # visualizziamo fs.19 rispetto a igD e igM: show(autoplot(fs.19[[1]],"PE-CF594-A")),show(autoplot(fs.19[[1]],"FITC-A"))
  # rispetto  a entrambi: show(autoplot(fs.19[[2]],"PE-CF594-A","FITC-A"))
  names(igm.gate)  <- names(igd.gate ) <- sampleNames(gs)
  
  IG <- fsApply(fs.19, function(x) {
    print(identifier(x))
    #igd.gate <- deGate (temp@flow.frame, channels.ind["IgD"],tinypeak.removal = 0.9,upper=F,alpha=0.9) 
    # codice per gate arbitrari:
    
    temp <- flowDensity(x, channels.ind[c("IgD","IgM")],position=c(T,NA),gates=c(igd.gate[identifier(x)],NA))
    igm.gate <- deGate ( temp@flow.frame, channels.ind["IgM"],percentile=NA,tinypeak.removal = .99,use.upper=T,upper=F,alpha=.8) 
    # if(igm.gate>=2.80){
    #   igm.gate<-igm.gate*0.70
    # }
    # igm.gate <- max(igm.gate,)
    quad.1 <- flowDensity(x, channels.ind[c("IgD","IgM")],position=c(F,F),gates=c(igd.gate[identifier(x)],igm.gate))
    quad.2 <- flowDensity(x,channels.ind[c("IgD","IgM")],position=c(T,F), gates=quad.1@gates)
    quad.3 <- flowDensity(x, channels.ind[c("IgD","IgM")],position=c(F,T), gates=quad.1@gates)
    quad.4 <- flowDensity(x, channels.ind[c("IgD","IgM")],position=c(T,T), gates=quad.1@gates)
    igm <- flowDensity(x, channels.ind[c("IgD","IgM")],position=c(NA,T), gates=quad.1@gates)
    return(list(quad1=quad.1,quad2=quad.2,quad3=quad.3, quad4=quad.4, igm=igm))
  })
  
  # ----- aggiorniamo il GatingSet
  
  IGigm_filter_list <- lapply(IG,function(x) x$igm@filter)
  IGigm.poly <- lapply(1:length(IGigm_filter_list),function(i){
    polygonGate(filterId = "IGigm",.gate=IGigm_filter_list[i])
  })
  names(IGigm.poly)<-sampleNames(gs)
  node_IGigm<-add(gs,IGigm.poly,parent="Bcells")
  recompute(gs)
  
  IGquad1_filter_list <- lapply(IG,function(x) x$quad1@filter)
  IGquad1.poly <- lapply(1:length(IGquad1_filter_list),function(i){
    polygonGate(filterId = "IGquad1",.gate=IGquad1_filter_list[i])
  })
  names(IGquad1.poly)<-sampleNames(gs)
  node_IGquad1<-add(gs,IGquad1.poly,parent="Bcells")
  recompute(gs)
  
  IGquad2_filter_list <- lapply(IG,function(x) x$quad2@filter)
  IGquad2.poly <- lapply(1:length(IGquad2_filter_list),function(i){
    polygonGate(filterId = "IGquad2",.gate=IGquad2_filter_list[i])
  })
  names(IGquad2.poly)<-sampleNames(gs)
  node_IGquad2<-add(gs,IGquad2.poly,parent="Bcells")
  recompute(gs)
  
  IGquad3_filter_list <- lapply(IG,function(x) x$quad3@filter)
  IGquad3.poly <- lapply(1:length(IGquad3_filter_list),function(i){
    polygonGate(filterId = "IGquad3",.gate=IGquad3_filter_list[i])
  })
  names(IGquad3.poly)<-sampleNames(gs)
  node_IGquad3<-add(gs,IGquad3.poly,parent="Bcells")
  recompute(gs)
  
  IGquad4_filter_list <- lapply(IG,function(x) x$quad4@filter)
  IGquad4.poly <- lapply(1:length(IGquad4_filter_list),function(i){
    polygonGate(filterId = "IGquad4",.gate=IGquad4_filter_list[i])
  })
  names(IGquad4.poly)<-sampleNames(gs)
  node_IGquad4<-add(gs,IGquad4.poly,parent="Bcells")
  recompute(gs)
  return(list(gs=gs,IG=IG))
}


# vediamo differenze di visualizzazzione
# plotDens(fs.19[[2]],channels=channels.ind[c("IgD","IgM")],main=identifier(x))#,nbin=3000
# abline(v=igd.gate[identifier(x)],h=igm.gate)
# 
# plotDens(getData(gs,"Bcells")[[1]],channels=channels.ind[c("IgD","IgM")],main=identifier(x))#,nbin=3000
# abline(v=igd.gate[2],h=igm.gate)




# ------------------------Gating from CD20+ to CD10- -----------------------
# according to the gate strategy il gate di cd10 e' sotto il 3
# show(autoplot(fs.20,"Alexa Fluor 700-A"))
# da tinypeak.removal 0.99 a 1/70,in questo modo escludo i picchi piccolissimi(piu' piccolo e' tinypeak.removal,piu' piccolo e' il picco che escludiamo)
#show(autoplot(fs.20[[1]],"Alexa Fluor 700-A"))
gating_cd20_to_cd10neg<-function(){
  cd10.gate <-averageGates(fsApply(fs.20,deGate,channels.ind[13],use.upper=T,upper=T,tinypeak.removal=0.9),2)
  names(cd10.gate)<-sampleNames(gs)
  CD10.27<- fsApply(fs.20, function(x){
    # temp1<- flowDensity(x,channels.ind[c(5,7)],position = c(NA,T),gates=c(NA,gd.gate[identifier(x)]))
    # cd3.gate.lo <- deGate (getflowFrame(temp1),channels.ind[7],upper=F,use.upper=T,tinypeak.removal = .9,alpha=.1)
    cd10.neg<-flowDensity(x,channels.ind[c(13,12)],position = c(F,NA),gates=c(cd10.gate[identifier(x)],NA))
    cd10.pos<-flowDensity(x,channels.ind[c(13,12)],position = c(T,NA),gates=c(cd10.gate[identifier(x)],NA))
    
    return(list(pos=cd10.pos, neg=cd10.neg))
  })
  
  # ----- aggiorniamo il GatingSet
  
  CD10pos.27_filter_list <- lapply(CD10.27,function(x) x$pos@filter)
  CD10pos.27.poly <- lapply(1:length(CD10pos.27_filter_list),function(i){
    polygonGate(filterId = "CD10pos.27",.gate=CD10pos.27_filter_list[i])
  })
  names(CD10pos.27.poly)<-sampleNames(gs)
  node_CD10pos.27<-add(gs,CD10pos.27.poly,parent="CD19.20pos")
  recompute(gs)
  
  CD10neg.27_filter_list <- lapply(CD10.27,function(x) x$neg@filter)
  CD10neg.27.poly <- lapply(1:length(CD10neg.27_filter_list),function(i){
    polygonGate(filterId = "CD10neg.27",.gate=CD10neg.27_filter_list[i])
  })
  names(CD10neg.27.poly)<-sampleNames(gs)
  node_CD10neg.27<-add(gs,CD10neg.27.poly,parent="CD19.20pos")
  recompute(gs)
  
  fs.10.neg <- getData(gs,"CD10neg.27")
  fs.10.pos <- getData(gs,"CD10pos.27")
  return(list(gs=gs,fs.10.neg=fs.10.neg,fs.10.pos=fs.10.pos,CD10.27=CD10.27))
}


# ------------------------ Gating from CD10- to Naive(CD27-,igd+) -----------------------
gating_cd10n_to_naive<-function(){
  # grafico rispetto a ogni dimensione: show(autoplot(fs.10.neg[1],"PE-Cy7-A")),show(autoplot(fs.10.neg[3],"PE-CF594-A"))
  # entrambe le dim: show(autoplot(fs.10.neg[[1]],"PE-Cy7-A","PE-CF594-A"))
  
  igd.gate <- averageGates(as.numeric(fsApply(fs.10.neg,function(x) return(deGate(x,channels.ind["IgD"],after.peak = F,percentile=NA,upper=F)[1]))),sd.coeff = 2)
  cd27.gate <- averageGates(fsApply(fs.ngran,deGate,channels.ind["CD27"]),sd.coeff=2)
  
  names(cd27.gate)  <- names(igd.gate) <-sampleNames(gs)
  
  naive <- fsApply(fs.10.neg, function(x){
    print(identifier(x))
    quad.1 <- flowDensity(x, channels.ind[c("CD27","IgD")],position=c(F,F),gates=c(cd27.gate[identifier(x)]*1.05,igd.gate[identifier(x)]))
    quad.2 <- flowDensity(x,channels.ind[c("CD27","IgD")],position=c(T,F), gates=quad.1@gates)
    quad.3 <- flowDensity(x, channels.ind[c("CD27","IgD")],position=c(F,T), gates=quad.1@gates)
    quad.4 <- flowDensity(x, channels.ind[c("CD27","IgD")],position=c(T,T), gates=quad.1@gates)
    return(list(quad1=quad.1,quad2=quad.2,quad3=quad.3, quad4=quad.4))
  })
  
  # ----- aggiorniamo il GatingSet
  
  filter_list <- lapply(naive,function(x) x$quad3@filter)
  poly <- lapply(1:length(filter_list),function(i){
    polygonGate(filterId = "naive",.gate=filter_list[i])
  })
  names(poly)<-sampleNames(gs)
  node_naive<-add(gs,poly,parent="CD10neg.27")
  recompute(gs)
  return(list(gs=gs,naive=naive,cd27.gate=cd27.gate,igd.gate=igd.gate))
}


# ------------------------ Gating from IGM+IGD+CD27- to Transitional cell-----------------------
gating_IGM_IGD_CD27_to_tran<-function(){
  transitional<- fsApply(fs.10.pos, function(x){
    f<- IG[[identifier(x)]]$igm
    cd27.igd <- flowDensity(f,  channels.ind[c("CD27","IgD")],position=c(F,T),gates=c(cd27.gate[identifier(x)],igd.gate[identifier(x)]))
    
    gate.1 <- deGate(cd27.igd, channels.ind["CD10"],twin.factor = .9,tinypeak.removal =.9,upper=T, alpha=.5, magnitude=.1)
    gate.2 <- deGate(cd27.igd, channels.ind["CD38"],twin.factor = .9,tinypeak.removal =.9,upper=F, alpha=.9)
    trans <- flowDensity(cd27.igd,  channels.ind[c("CD10","CD38")],position=c(T,T),gates=c(gate.1,gate.2))
    # lines(dr.neg@filter,type="l",lwd=2)
    return(list(trans=trans,parent=cd27.igd))
  })
  
  # ----- aggiorniamo il GatingSet
  filter_list <- lapply(transitional,function(x) x$trans@filter)
  poly <- lapply(1:length(filter_list),function(i){
    polygonGate(filterId = "trans",.gate=filter_list[i])
  })
  names(poly)<-sampleNames(gs)
  node_trans<-add(gs,poly,parent="CD10neg.27")
  recompute(gs)
  return(list(gs=gs,transitional=transitional))
} 

# ------------------------ Gating  from live to blast cells-----------------------
# proportion manual count: 0.005496685,0.003165626,0.1336438
# show(autoplot(fs.live[1],"PE-A","SSC-A"))
# show(autoplot(fs.live,"SSC-A"))
# show(autoplot(fs.live[1],"PE-A"))
gating_live_to_blast<-function(){
  ss.ind <- grep("SSC-A",colnames(fs[[1]]))
  Blast<- fsApply(fs.live, function(x){
    gate.1 <- deGate(x,ss.ind,upper=F)
    temp <- flowDensity(x,c(channels.ind[4],ss.ind),position = c(NA,F),gates=c(NA,gate.1))
    gate.2 <- deGate(temp,channels.ind[4],upper=T,use.upper=T,alpha = .02)
    blast <-flowDensity(x,c(channels.ind[4],ss.ind),position = c(T,F),gates=c(gate.2,gate.1))
    #print(gate.2)
    return(blast)
  })
  
  # ----- aggiorniamo il GatingSet
  
  Blast_filter_list <- lapply(Blast,function(x) x@filter)
  Blast.poly <- lapply(1:length(Blast_filter_list),function(i){
    polygonGate(filterId = "Blast",.gate=Blast_filter_list[i])
  })
  names(Blast.poly)<-sampleNames(gs)
  node_Blast<-add(gs,Blast.poly,parent="Live")
  recompute(gs)
  return(list(gs=gs,Blast=Blast))
}




#--------  Function to import and set up the dataset of manual counts

pre_prop_manual_data <- function(){
  manual_data<- read.csv(file="/home/rstudio/data/PNG B cell panel (Extended)/PNG_Pilot_B_cell_counts_export_(all populations).csv",header=TRUE,sep=",")

  #------- riassegnamo nomi colonne del dataset manuale
  pop_names_manual <- colnames(manual_data)
  pop_names_manual[which(pop_names_manual=="X")]<- "Sample"
  pop_names_manual[which(pop_names_manual=="Blasts")]<- "Blast"
  pop_names_manual[which(pop_names_manual=="CD66.neg.CD14.neg")]<- "NonGran"
  pop_names_manual[which(pop_names_manual=="CD19...B.cells")] <- "Bcells"
  pop_names_manual[which(pop_names_manual=="CD19..CD20.")] <- "CD19.20pos"
  pop_names_manual[which(pop_names_manual=="CD27.CD10.")] <- "CD10pos.27"
  pop_names_manual[which(pop_names_manual=="CD27.CD10..")] <- "CD10neg.27"
  pop_names_manual[which(pop_names_manual=="CD27.CD10..CD27..IgD.")] <- "naive"
  pop_names_manual[which(pop_names_manual=="Plasmablats")] <- "PB"
  pop_names_manual[which(pop_names_manual=="CD19.positive.gM..IgD.")] <- "IGquad3"
  pop_names_manual[which(pop_names_manual=="CD19.positive.IgM..IgD.")] <- "IGquad4"
  pop_names_manual[which(pop_names_manual=="CD19.positive.gM...IgD.")] <- "IGquad2"
  pop_names_manual[which(pop_names_manual=="CD19.positive.gM..IgD..")] <- "IGquad1"
  colnames(manual_data) <- pop_names_manual
  #------ adesso convertiamo le conte in percentuali delle colonne delle popolazioni
  # manual_data_perc <- manual_data[,2:length(colnames(manual_data))]/manual_data[,"Live cells"] # da sei se vogliamo includere le pop solo a partire dalle live,da 1 se tutte
  #--- creiamo dataframe che associa a ogni nome di una pop la sua parente
  list_child_vs_parent<-lapply(1:length(getNodes(gs)), function(i){
    splitted_string <-str_split(getNodes(gs)[i],"/")
    child_pop <- splitted_string[[1]][length(splitted_string[[1]])]
    parent_pop <- splitted_string[[1]][length(splitted_string[[1]])-1]
    return(list(child_pop=child_pop,parent_pop=parent_pop))
  })

  list_all_pop_freq <- lapply(1:length(list_child_vs_parent), function(i){
    child_pop_x <- list_child_vs_parent[[i]]$child_pop

    if(child_pop_x %in% colnames(manual_data)){
      parent_pop_x <- list_child_vs_parent[[i]]$parent_pop
      if(parent_pop_x %in% colnames(manual_data)){
        pop_x_freq <- manual_data[,child_pop_x]/manual_data[,parent_pop_x]

        return(list(pop_x_freq=pop_x_freq,name=child_pop_x))
      }
    }
  })

  indx<-which(list_all_pop_freq=="NULL")
  list_all_pop_freq<-list_all_pop_freq[-indx]

  names<-sapply(1:length(list_all_pop_freq),function(i){
    return(list_all_pop_freq[[i]]$name)
  })
  list_all_pop_freq_only<-lapply(1:length(list_all_pop_freq),function(i){
    return(list(pop=list_all_pop_freq[[i]]$pop_x_freq))
  })
  df <- as.data.frame(list_all_pop_freq_only)
  colnames(df)<-names
  manual_data_freq_on_parents <- df
  samples<-manual_data$Sample
  samples_split<-strsplit(as.character(samples),": ")
  samples<-sapply(1:length(samples_split),function(i){
    return(samples_split[[i]][2])
  })
  manual_data_freq_on_parents<-cbind(samples,manual_data_freq_on_parents)
  return(manual_data_freq_on_parents)
}





####################################### Execute functions ##################################
main<-function(n_sample="All",path_comp_matrix,path_dir_files,path_fcs_files_fixing=NULL,plots=T,flowjo_wsp=F){
  #Eseguiamo le funzioni create
  # resettiamo il gs
  #gs<- get_GatingSet()
  #----- preprocessing
  # importiamo funzioni utili al pre-processing
  setwd("/home/rstudio/Code_Bcells_Myeloid_cells")
  path_utils <<- paste0("./Utils.R")
  source(path_utils)
  import_flowPrep()
  output_import_files<<-import_files_v2(n_sample=n_sample,files_path=path_dir_files)
  
  fs<<-output_import_files$fs
  FCSfiles_path_list <<- output_import_files$FCSfiles_path_list
  path.output <<- create_path_output_dir("Bcells_PNG")
  channels.ind<<-find_markers_indices_bcells()
  gs<<-get_GatingSet()
  # ---- gating dei marginal events
  output<<-margin_gating()
  fs.marg<<-output$fs.marg
  gs<<-output$gs
  # ----- compensation e trasformation
  comp<<-import_comp_matrix_v2(path_comp_matrix=path_comp_matrix)
  out<<-comp_and_transform()
  gs<<-out$gs
  fs.marg<<-out$fs.marg
  #----- cleaning
  output_prova <<- cleaning()
  fs.clean <<- output_prova$fs.clean
  gs <<- output_prova$gs
  Clean <<- output_prova$Clean
  #--- gating time to singlets
  output_prova_2 <<- gating_time_to_singlets()
  gs <<- output_prova_2$gs
  fs.sngl <<- output_prova_2$fs.sngl
  Singlet <<- output_prova_2$Singlet
  # Rm("Singlets",gs)
  #plotGate(gs,"Singlets")
  #------- Size e beads
  output_prova_3 <<- gating_singlets_to_size_beads()
  fs.size <<- output_prova_3$fs.size
  gs <<- output_prova_3$gs
  size <<- output_prova_3$size
  beads <<- output_prova_3$beads
  
  # fs.size  <<- fsApply(fs.size, function(x) return(compensate(x, as.matrix(comp))))
  # fs.size <<-transform.flow(fs.size,remove.outliers = T,sd.coeff = 3.5,trans.chans = c(7:19))
  # show(autoplot(fs.size[[1]],"APC-eF780-A","SSC-A"))
  # Rm("Size",gs)
  # Rm("Beads",gs)
  # plotGate(gs,"Size")
  # plotGate(gs,"Beads")
  #---- Live
  output_prova_4 <<- gating_size_to_live()
  fs.live<<-output_prova_4$fs.live
  gs<<-output_prova_4$gs
  live<<-output_prova_4$live
  # Rm("Live",gs)
  # plotGate(gs,"Live")
  #------ to Granulocytes e non gran
  output_prova_5 <<- gating_live_to_non_Granulocytes()
  gs <<- output_prova_5$gs
  fs.ngran <<- output_prova_5$fs.ngran
  nonGran <<- output_prova_5$nonGran
  #plotGate(gs,"NonGran")
  # Rm("NonGran",gs)
  #--- Bcells(CD19+)
  output_prova_6 <<-gating_nonGran_to_Bcells()
  gs<<-output_prova_6$gs
  fs.19<<-output_prova_6$fs.19
  Bcell<<-output_prova_6$Bcell
  # Rm("Bcells",gs)
  #plotGate(gs,"Bcells")
  #--- bcells to PC
  output_prova_7 <<- gating_bcells_to_pc()
  gs<<-output_prova_7$gs
  PC<<-output_prova_7$PC
  # Rm("PC",gs)
  #plotGate(gs,"PC")
  #--- bcells to PB
  output_prova_8<<-gating_bcells_to_PB()
  gs <<- output_prova_8$gs
  PB <<- output_prova_8$PB
  #Rm("PB",gs)
  #plotGate(gs,"PB")
  # --- Bcells to CD19/20
  output_prova_9<<-gating_bcells_to_CD20p()
  gs<<-output_prova_9$gs
  fs.20<<-output_prova_9$fs.20
  CD19.20<<-output_prova_9$CD19.20
  # plotGate(gs,"CD19.20pos")
  # plotGate(gs,"CD19.20neg")
  # Rm("CD19.20pos",gs)
  # Rm("CD19.20neg",gs)
  # --- Bcells to igd/igm
  output_prova_10<<-gating_bcells_to_igm_igd()
  gs<<-output_prova_10$gs
  IG<<-output_prova_10$IG
  # plotGate(gs,"IGquad1")
  # plotGate(gs,"IGquad2")
  # plotGate(gs,"IGquad3")
  # plotGate(gs,"IGquad4")
  # Rm("IGigm",gs)
  # Rm("IGquad1",gs)
  # Rm("IGquad2",gs)
  # Rm("IGquad3",gs)
  # Rm("IGquad4",gs)
  #------- cd20 to cd10neg
  output_prova_11<<-gating_cd20_to_cd10neg()
  gs<<-output_prova_11$gs
  fs.10.neg<<-output_prova_11$fs.10.neg
  fs.10.pos<<-output_prova_11$fs.10.pos
  CD10.27<<-output_prova_11$CD10.27
  # plotGate(gs,"CD10pos.27")
  # plotGate(gs,"CD10neg.27")
  # Rm("CD10pos.27",gs)
  # Rm("CD10neg.27",gs)
  #----- cd10neg to Naive
  output_prova_12<<-gating_cd10n_to_naive()
  gs<<-output_prova_12$gs
  naive<<-output_prova_12$naive
  cd27.gate<<-output_prova_12$cd27.gate
  igd.gate<<-output_prova_12$igd.gate
  #plotGate(gs,"naive")
  #Rm("naive",gs)
  # ---- to trans
  output_prova_13<<-gating_IGM_IGD_CD27_to_tran()
  gs<<-output_prova_13$gs
  transitional<<-output_prova_13$transitional
  #plotGate(gs,"trans")
  #Rm("trans",gs)
  #----- from live to blast
  output_prova_14<<-gating_live_to_blast()
  gs<<-output_prova_14$gs
  Blast<<-output_prova_14$Blast
  #Rm("Blast",gs)
  #plotGate(gs,"Blast")
  #---- Generate final Data
  if(flowjo_wsp==T){
    generate_final_data()
  }
  #---- Generates plots of each Gating step for each sample
  if(plots==T){
    generate_plots_Bcells()
  }
  return(gs)
}





#------------- section not usefull in fina version: comparison with manual counts --------------------
#---- generate correlation plot
# manual_data_freq_on_parents<-pre_prop_manual_data()
# tab_stat_counts_freq<-pre_prop_automated_data()
# final_df<-make_final_dataframe()
# # fix: IGquad gates(gate near the peal of the density),CD19.29pos(gate orizzontale piu' in alto),
# # Cd10neg.27,gate verticale piu' a sinistra
# corr_plot<-make_corr_plot()
# corr_plot<-make_corr_plot(remove_p = "IGquad4|IGquad3|IGquad2|IGquad1|CD19.20neg|CD19.20pos|CD10pos.27|CD10neg.27")
# 
# show(corr_plot)
# # ----- generate boxplots
# boxplots <- generate_box_plots()
main("All",path_comp_matrix ="/home/rstudio/data/Sample info/Bcell-comp.csv", 
     path_dir_files = "/home/rstudio/data/PNG B cell panel (Extended)/PNG B cells expended panel 01032018/")

