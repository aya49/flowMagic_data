

set.seed(40)

#---------------- create margin gates -------------------------------
margin_gating <- function(){
  gs <-get_GatingSet()
  # otteniamo il flowSet dal GatingSet con getData(). Ricorda che il flowSet e' costituito da un solo FlowFrame in questo caso.
  fs <- getData(gs)
  # Usiamo la funzione removeMarging del file flowPrep.R sui channels FCS-A,SSC-A che sono i channels che servono per identificare gli eventi marginali
  # L'output sono i gates,ovvero il loro threshold
  margin.gates <- fsApply(fs,removeMargins,c("FSC-A","SSC-A"),return.gate=T)
  # margin.gates e' formato da due colonne,la prima indica il threshold per FCS-A,la seconda indica il Threshold per SSC-A
  # usiamo la funzione lapply avendo come input il numero di righe di margin.Gates (ovvero solo una riga) 
  # e settiamo le coordinate di un rettangleGate in base ai threshold per il gating dei margins(il boundary lo facciamo partire da -1)
  # Nota: fino al threshold indicato da margin.gates sono non margin cells,dopo il threshoold indicato da margin.gates abbiamo le margin cells.
  rg <- lapply(1:nrow(margin.gates), function(x) return(rectangleGate(filterId="Margin", "FSC-A"=c(-1, margin.gates[x,1]),
                                                                      "SSC-A"=c(-1, margin.gates[x,2]))))
  # assegnamo come nome di rg(che attualmente e' NULL) il nome dei samples di fs.
  names(rg) <-sampleNames(fs)
  # Aggiugiamo questo popolazione/nodo/gate al GatingSet(e quindi al gating tree dei suoi GatingHierarchy,in questo caso solo 1)
  nodeID1<-add(gs, rg)#,"Margin")
  # add presenta la seguente struttura:
  # add(wf, action, ...) dove:
  # wf = A GatingHierrarchy or GatingSet
  # action = A filter or a list of filters to be added to the GatingHierarchy or GatingSet.
  # Eseguiamo effettivamente il Gating sul flowSet del GatingSet in base alle informazioni che abbiamo appena caricato.
  recompute(gs)
  # prendiamo il flowSet del Gating in corrispondenza del nodo Margin. 
  # Cioe' in pratica prendiamo le cell ripulite dai margins events,quindi la subpopolazione appena creata.
  fs_margin<-getData(gs,"Margin")
  return(list(fs_margin=fs_margin,gs=gs))
}
#-------------- Compensation Trasformation----------------------




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
  fs.marg<-getData(gs,"Margin")
  return(list(gs=gs,fs.marg=fs.marg))
}
#------------------ cleaning ---------------------------
cleaning<-function(){
  clean.inds <- fsApply(fs.marg, function(f){
    # f1 <- compensate(f,comp)
    # f1 <- estimate.logicle(fs.raw = as(f1,"flowSet"),med=F,estimate=T,trans.chans = c(7:19))
    #f1 <- f1[[1]]
    res.tvsf <- flowCut(f=f, Segment=10000, Channels=NULL, Directory= paste0(path.output,"/Cleaning/"), 
                        FileID= paste0(identifier(f),"_"), Plot="All")
    return(res.tvsf)
  })

  clean.clust <- lapply(1:length(clean.inds), function(x){
    # ripetiamo il numero 0 un numero di volte pari alle righe degli eventi presenti nel flowFrame del flowSet dei non marginal events.
    vec<-rep(0,nrow(fs.marg[[x]]))
    # la lunghezza dello slot degli indici dell'output di flowCut è maggiore di 0. Quindi significa che sono stati trovati eventi da "rimuovere".
    if (length(clean.inds[[x]]$ind)>0)
    {
      # sostituiamo tutti elementi del vettore in corrispondenza degli indici dello slot ind di clean.inds con 1.
      # questo significa che gli elementi del vettore che non rientrano negli indici di clean.inds rimangono 0 e sono le "good cells",
      # perche'è [clean.inds[[x]]$ind riporta gli indici delle cellule da rimuovere ovvero quelle "cattive"
      vec[clean.inds[[x]]$ind]<-1
    }
    # assegnamo alla classe factor gli elementi di vec
    vec <- as.factor(vec)
    # come livelli del fattore vec assegnamo 0 e 1. 0 = good cells, 1= bad cells
    levels(vec) <- c("0", "1")
    return(vec)
  })
  # dunque clean.clust e' un vettore fattorizzato che contiene due livelli indicanti se l'evento in quella posizione e' data tenere(clean events,quindi fattore 0)
  #oppure da non considerare(bad event,quindi fattore = 1)
  #assegnamo sempre al nome dell'oggetto clean.clust il nome del sample del flowSet.
  names(clean.clust)<-sampleNames(gs)
  # aggiungiamo il  nodo delle cellule pulite al gating tree del GatingSet.
  # Nota che stavolta abbiamo aggiunto un vettore fattorizzato(praticamente un vettore di clusters) che agisce come filtro fra gli gli eventi della popolazione madre(margin)
  # e la nuova popolazione(Clean). Prima,invece nel caso del margin node,avevamo aggiunto un oggetto gate filter con gli estremi dei valori degli eventi della popolazione figlia.
  #  ATTENZIONE: aggiungendo un vettore di fattorizzazione con due livelli,creiamo due nodi che prenderanno il nome "nome_nodo_1" e "nome_nodo_2".
  # Quindi generiamo in un colpo solo due popolazioni,quella solo con eventi puliti(Clean_0) e quella solo con i cattivi eventi(Clean_1)
  add(gs,clean.clust,parent = "/Margin", name = "Clean")
  recompute(gs)

  fs.clean <- getData(gs, "Clean_0")
  return(list(clean.inds=clean.inds,fs.clean=fs.clean,gs=gs))
}

# ------------------ Gating Time to get Singlets -------------------------------------------
gating_Time_to_singlets<-function(){
  singlets_poly_list <- fsApply(fs.clean, function(f) {
    print(identifier(f))
    channels <- c("FSC-A", "FSC-H")
    rot <- rotate.data(f, channels,theta = -pi/4)$data
    gate.2 <- deGate(rot, channels[1],percentile=NA,tinypeak.removal = .98,twin.factor = .2,magnitude = .05, upper=T, alpha=.007,verbose = F)*1.06
    singlets<- rotate.fd(flowDensity(rot,channels,position = c(F,NA),gates=c(gate.2,NA),verbose = F),angle = -pi/4)
    singlets<-singlets@filter
    sngl.poly <- polygonGate(filterId = "Singlets",.gate=singlets)
    return(list(sngl.poly=sngl.poly,singlets=singlets))
  })
  list_singlets_poly <- lapply(1:length(gs),function(i){
    singlets_poly_list[[i]]$sngl.poly
  })
  list_singlets <- lapply(1:length(gs),function(i){
    singlets_poly_list[[i]]$singlets
  }) # lista che contiene i filters del population object dei singlets per ogni sample
  #assegnamo un nome a questo oggetto pari al nome del sample
  names(list_singlets_poly)<-sampleNames(gs)
  # aggiungiamo il filtro al gating tree
  nodeID0<-add(gs,list_singlets_poly,parent="Clean_0")
  recompute(gs)
  # estraiamo popolazione dei Singlets
  fs.sngl<-getData(gs,"Singlets")
  return(list(fs.sngl=fs.sngl,gs=gs,list_singlets=list_singlets))
}


# --------------------------------- Gating singlets to get Beads, Size -------------------------------------

gating_singlets_to_Size_beads <- function(){

  beads_size_poly_list <- fsApply(fs.sngl, function(x){
    #----- prima facciamo il gating dei size
    print(paste("Gating size of ",identifier(x)))
    
    #x<-removeMargins(x,c("FSC-A","SSC-A"),neg = 500,sens = .999)
    x.rot <- rotate.data(data = x,chans =c("FSC-A","SSC-A"),theta=-pi/2.5)$data
    peaks <-getPeaks(x.rot, "SSC-A",tinypeak.removal = .02)$Peaks
    if(peaks[which.min(peaks-55000)]>70000)
      ss.lo <-  deGate(x.rot, "SSC-A", upper=F,use.upper=T, percentile=NA,alpha=.1)
    else
      ss.lo <- deGate(x.rot, "SSC-A", upper=NA, percentile=NA,alpha=0.025,all.cut=T,tinypeak.removal = .02,n.sd=.4)[1]*.95
    
    temp <-flowDensity(x.rot,channels =c("FSC-A","SSC-A"), position=c(NA,T),gates=c(NA,ss.lo))
    temp@flow.frame <- rotate.data(getflowFrame(temp),chans =c("FSC-A","SSC-A"),theta =pi/2.5)$data
    ss.hi <- deGate(getflowFrame(temp), "SSC-A", upper=T, percentile=NA,alpha=0.1,all.cut=T,tinypeak.removal = 1/40)
    ss.hi <-ss.hi[which.min(abs(ss.hi-170000))]
    size <- flowDensity(temp,channels =c("FSC-A","SSC-A"), position=c(NA,F),gates=c(NA,ss.hi))
    size@proportion<-size@cell.count/nrow(x)*100
    size.poly <- polygonGate(filterId = "All cells",.gate=size@filter)
    #--- adesso facciamo il gating dei beads
    print(paste("Gating beads of ",identifier(x)))
    temp <-flowDensity(x,channels =c("FSC-A","SSC-A"), position=c(NA,T),gates=c(NA,size@gates[2]))
    gate <- deGate(temp, "FSC-A", percentile=NA,alpha=.01,use.upper=T,upper=T,tinypeak.removal = .8)
    temp2 <- flowDensity(temp,channels =c("FSC-A","SSC-A"), position=c(F,NA), gates=c(gate,NA))
    bead <- flowDensity(temp2,channels =c("FSC-A","SSC-A"), position=c(F,NA), use.percentile=c(T,F),percentile=c(.999,NA),ellip.gate=T)
    bead@proportion<-bead@cell.count/nrow(x)*100

    beads.poly <- polygonGate(filterId = "Beads",.gate=bead@filter)
    return(list(size.poly=size.poly,size=size,bead=bead,beads.poly=beads.poly))
  })


  
  beads_poly_list<-lapply(1:length(gs),function(i){
    beads_size_poly_list[[i]]$beads.poly
  })
  size_poly_list<-lapply(1:length(gs),function(i){
    beads_size_poly_list[[i]]$size.poly
  })
  
  list_beads<-lapply(1:length(gs),function(i){
    beads_size_poly_list[[i]]$bead
  })

  # assegnamo il nome
  names(size_poly_list)<-sampleNames(gs)
  # aggiungiamo il nodo al gating tree
  nodeID0<-add(gs,size_poly_list,parent="Singlets")
  recompute(gs)
  # assegnamo il nome di questa lista al nome del sample
  names(beads_poly_list)<-sampleNames(gs)
  # aggiungiamo il nodo beads al gating tree
  nodeID0<-add(gs,beads_poly_list,parent="Singlets")
  recompute(gs)
  # otteniamo il flowSet della size population
  fs.size<-getData(gs,"All cells")
 
  return(list(fs.size=fs.size,gs=gs,beads_poly_list=beads_poly_list,size_poly_list=size_poly_list,list_beads=list_beads))
}
  



# ------------------ gating delle live a partire dalla Size population(All cells population) ---------

gating_size_to_live<- function(){

  live.gate <- averageGates(fsApply(fs.size,function(x) tail(deGate(x,"APC-eF780-A",upper=T,use.upper = T, tinypeak.removal = .5,alpha=0.05),1)),sd.coeff = 3.5)
  names(live.gate) <- sampleNames(fs.size)
  live_poly_list <- fsApply(fs.size, function(x){
    print(paste("Gating beads of ",identifier(x)))
    live <-flowDensity(x,channels =c("APC-eF780-A","SSC-A"), position=c(F,NA),gates=c(live.gate[identifier(x)],NA))
    # live.up <- deGate(temp,"APC-eF780-A",upper=F)
    #live <-flowDensity(x,channels =c("APC-eF780-A","SSC-A"), position=c(F,NA),gates=c(live.up,NA))
    live.poly <- polygonGate(filterId = "Live cells",.gate  =live@filter)
    return(live.poly)
  })
  # creiamo il poligon base in base al filter object di flowDensity
  names(live_poly_list)<-sampleNames(gs)
  # aggiungiamo il nodo delle live cells a partire dal nodo "All cells"
  nodeID0<-add(gs,live_poly_list,parent="All cells")
  recompute(gs)
  fs.live <- getData(gs, "Live cells")
  return(list(fs.live=fs.live,gs=gs,live_poly_list=live_poly_list))
}


# ---------------------- Gating Live to get Granulocytes e CD45+CD66- -------------------------------------------
gating_live_to_granulocytes<-function(){
  # adesso performiamo il gating dalle live per ottenere i Granulocytes(ovvero la pop CD66+) e la pop CD45+CD66-
  # ---------- Troviamo prima la pop dei granulociti 
  # adjust.dens :The smoothness of density in [0,Inf] to be used in density(.). The default value is 1
  # all.cuts restiuisce tutt i possibili thresholds
  # noi però prendiamo solo il primo threshold infatti scriviamo alla fine: [1]
  # visualizziamo la pop fs.live rispetto a CD66: show(autoplot(fs.live[[1]],"BV711-A","V450-A"))
  list_poly_gran_cd45<-lapply(1:length(gs),function(i){
    CD66.gate <- deGate(fs.live[[i]], channel = c(fluorochrome.chans['CD66']), all.cuts = T, adjust.dens = .8)[1]
    # e' pari a 2.577935 nel file HIPC PNG_Tube_001.fcs ,quindi siamo ancora un po' lontani dal 3.
    # in realtà all.cuts restituisce un solo tresholds che corrisponde al threshold calcolato con le impostazioni di default,ovvero percentile = NA
    # se questo e' maggiore del treshold calcolato usando l'85th percentile allora significa che e' un treshold che praticamente sta dopo due picchi
    # e non in mezzo. Allora in questo caso calcoliamo il threshold rimuovendo  i picchi più piccoli, 
    # in questo caso rimuoviamo i picchi che sono pari a 1/50 della densità e se nonostante questo ci sono ancora più di due picchi
    # troviamo un treshold che cmq divida la popolazione in due parti,infatti settiamo bimodal = True
    print(identifier(fs.live[[i]]))
    if (CD66.gate> deGate(fs.live[[i]],  channel = c(fluorochrome.chans['CD66']),use.percentile=T,percentile=.85))
      CD66.gate <- deGate(fs.live[[i]], channel = c(fluorochrome.chans['CD66']), tinypeak.removal =1/50, bimodal=T)
    # bidmodal = Logical. If TRUE, it returns a cutoff that splits population closer to 50-50, when there are more than two peaks. 
    # calcoliamo i gate boundaries e una popolazione temporenea considerando dunque il treshold appena calcolato. Noi stiamo interessati alla popolazione CD66+,CD45 per il momento non ci interessa
    temp <- flowDensity(fs.live[[i]] , channels = c(fluorochrome.chans['CD66'], fluorochrome.chans['CD45']), position = c(T, NA), gates = c(CD66.gate, NA))
    # plotDens(fs.live[[15]],c(fluorochrome.chans['CD66'],fluorochrome.chans['CD45']))
    # Dunque prediamo questa popolazione CD66+ e calcoliamo il treshold nella coda della distribuzione di CD66 quindi upper=F
    CD66.gate.up<- deGate(temp,fluorochrome.chans['CD66'],upper=F,alpha = 0.8)*0.93
    # e' pari a 3.045182 nel file HIPC PNG_Tube_001.fcs ,quindi effettivamente ci stiamo spostati un po' piu' a destra.
    # guardiamo la densità di questa pop temporanea rispetto  a CD66. autoplot(temp@flow.frame,"BV711-A") 
    # Ricorda che autoplo accetta solo il nome dei channels o i lori indici e non i nomi dei markers
    # come vedi adesso c'è un enorme grande picco centrale, mentre la CD66 precedente aveva due picchi,noi in pratica vogliamo quello di destra per ora.
    # con il treshold precedente abbiamo effettivamnete preso la pop che comprende il picco di destra,ma vogliamo avvicinarci di più e infatti abbiamo ricalcolato un nuovo threshold
    # facendolo corrispondere a valle del picco.
    #  a questo punto creiamo i gate boundaries della popolazione CD66+ con questo nuovo treshold a valle,sempre partendo però dalla popolazione delle Live originale
    # ovvero quella che ancora i CD66 con i due picchi. Adesso infatti conosciamo il trheshold che mi porta più vicino al pop del picco di destra.
    # abbiamo usato la popolazione CD66+ precedente soltanto in maniera temporanea per calcolare un threshold più specifico per la popolazione dei granulocites.
    granulocytes<- flowDensity(fs.live[[i]], channels = c(fluorochrome.chans['CD66'], fluorochrome.chans['CD45']), position = c(T, NA), gates = c(CD66.gate.up, NA))
    # Abbiamo selezionato dunque la popolazione dei granulociti,ovvero le cellule che sono CD66+ (valori di CD66 intorno al picco di 3) in alta concentrazione (picco molto alto)
    # vogliamo adesso prendere i granulociti(CD66+) che sono anche CD45+
    # guardiamo come e' organizatta la densità rispetto a CD45 della pop CD66+: show(autoplot(granulocytes@flow.frame,fluorochrome.chans['CD45']))
    #  c'è un grosso picco a sinistra ( anche se non molto alto)
    # prendiamo la popolazione di questo picco (ricorda che siamo dentro la popolazione granulocytes ovvero quella con CD66+,picco di CD66 di destra molto alto)
    # prendiamo il treshold a destra del picco
    CD45granulo.upper.gate <- deGate(granulocytes, channel = c(fluorochrome.chans['CD45']), use.upper = T, upper = T, alpha = 0.05) + 0.05
    # prendiamo il treshold a sinistra del picco
    CD45granulo.lower.gate <- deGate(granulocytes, channel = c(fluorochrome.chans['CD45']), use.upper = T, upper = F)
    # prendiamo la popolazione che include solo il picco. Quindi la densità fra questi due treshold.
    # Quindi prendiamo i gate boundaries e la popolazione a destra del treshold più basso
    granulocytes<- flowDensity(granulocytes, channels = c(fluorochrome.chans['CD66'], fluorochrome.chans['CD45']), position = c(NA, T), gates = c(NA, CD45granulo.lower.gate))
    # da questa popolazione prendiamo la popolazione con densità di CD45 a sinistra del treshold più alto.
    granulocytes <- flowDensity(granulocytes, channels = c(fluorochrome.chans['CD66'], fluorochrome.chans['CD45']), position = c(NA, F), gates = c(NA, CD45granulo.upper.gate))
    # calcoliamo con la funzione density() la density sulla matrice di espressione del flowFrame di quest'ultima pop (del channel CD45)
    maxDens <- density(getflowFrame(granulocytes)@exprs[, c(fluorochrome.chans['CD45'])])
    # prendiamo i valori di espressione del flowFrame(l'input x di density) in corrispondenza del picco di densità calcolato da density()(massima densità y)
    granulocytes.CD45peak <- maxDens$x[which.max(maxDens$y)]
    # Dunque la popolazione che abbiamo ottenuto finora,e' una popolazione che ha CD66+ e CD45+ anche se il picco di CD45 presenta valori non altissimi(vicino all'1)
    # abbiamo salvato il valore della punta del picco della densità CD45 della pop CD66+
    
    #------- Adesso troviamo la pop CD45+CD66-
    # generiamo la popolazione temporanea CD66- usando il primo treshold. Infatti nella prima temp population 
    # avevamo preso la pop a destra del treshold,ora prendiamo quella a sinistra ( picco con valori più bassi di CD66)
    # visualizziamo pop live rispetto a CD45 e CD66: show(autoplot(fs.live[[3]],"BV711-A")),show(autoplot(fs.live[[3]],"V450-A"))
    temp <- flowDensity(fs.live[[i]], channels = c(fluorochrome.chans['CD66'], fluorochrome.chans['CD45']), position = c(F, NA), gates = c(CD66.gate, NA))
    # calcoliamo un altro treshold su questa popolazione temporanea CD66- nella testa della distribuzione
    temp <- flowDensity(fs.live[[i]], channels = c(fluorochrome.chans['CD66'], fluorochrome.chans['CD45']), position = c(F, NA), gates = c(CD66.gate, NA))
    # calcoliamo un altro treshold su questa popolazione temporanea CD66- nella testa della distribuzione
    CD66.gate.lo<- deGate(temp,fluorochrome.chans['CD66'],tinypeak.removal = 0.9,upper=T,alpha=0.08)*1.05
  
    # prendiamo la popolazione a sinistra di questo treshold,quindi in pratica la pop che comprende quasi interamente il picco di sinistra di CD66,dunque una pop CD66-
    CD66negCD45pos<-flowDensity(fs.live[[i]], channels = c(fluorochrome.chans['CD66'], fluorochrome.chans['CD45']), position = c(F, NA), gates = c(CD66.gate.lo, NA))
    # osserviamo la densità di CD45 per questa popolazione CD66-:show(autoplot(CD66negCD45pos@flow.frame,"V450-A"))
    # Nota che c'è un errore nel flow.frame il marker il channel V450-A indica CD45 e non CD16
    # Come vedi c'è una serie di picchi verso destra, quindi valori di CD45 tendenzialmente molto più alti dei CD45 della pop CD66+ che abbiamo visto prima.
    # Su questa popolazione CD66- calcoliamo tutti i possibili thresholds per la densità di CD45 ( come vedi ci sono molti picchi quindi molti possibili thresholds)
    CD45.gate <- deGate(CD66negCD45pos, channel = c(fluorochrome.chans['CD45']), use.upper=F,upper=F, alpha=0.6)*0.85
    #CD45.gate <- max(CD45.gate[which(CD45.gate < granulocytes.CD45peak)])
    
    
    # fra questi treshold selezioniamo quelli che sono inferiori al valore esatto del picco(la punta proprio) 
    # di densità di CD45 per la popolazione CD66+ che abbiamo trovato prima. 
    # Prendiamo il treshold massimo fra questi,quindi quello più a destra.
    # ti ricordo che il picco di Cd45 della pop CD66+ era molto a sinistra,intorno a 1, rispetto a questi picchi di Cd45 della pop Cd66-
    # e noi abbiamo preso il valore più a destra di tutti i valori di sinistra di questo picco
    
    # tuttavia se il treshold selezionato e' maggiore 2.2 e il primo dei picchi della popolazione CD66- 
    # per il channels CD45 è maggiore di 1.1....
    # getPeaks = Find all peaks in density along with their indices
    # Allora e' un treshold troppo spostato  destra,allora riassegnamo questo threshold come treshold più a sinistra (upper=F) di tutti i picchi.
    # if (CD45.gate >2.2 & getPeaks(CD66negCD45pos, fluorochrome.chans['CD45'],tinypeak.removal = 0.001)$Peaks[1]>1.1)
    #   CD45.gate <- deGate(CD66negCD45pos, channel = c(fluorochrome.chans['CD45']),upper=F, use.upper = T, alpha=.1)
    # ritorniamo dunque alla popolazione live iniziale. Stavolta conosciamo i treshold ideali per avere CD66- e CD45+.
    # come threshold per CD66- scegliamo il treshold più preciso(cioè estratto dalla pop temporanea) che avevamo usato prima per calcolare la pop CD66-.
    # come treshold per CD45+ usiamo il treshold ideale appena trovato sulla pop CD66 che come abbiamo già detto e' stato calcolato su 
    # una pop che presenta Cd45 tendenzialmente più alti rispetto a ciò che abbiamo visto per la popolazione CD66+
    CD66negCD45pos<- flowDensity(fs.live[[i]], channels = c(fluorochrome.chans['CD66'], fluorochrome.chans['CD45']), position = c(F, T), gates = c(CD66.gate.lo, CD45.gate))
    # Perfetto,abbiamo dunque i nostri Gate boundaries con la popolazione CD66-CD45+ che volevamo.
    # adesso costruiamo i Gate objects(sfruttando i boundaries appena creati) dei granulociti e di CD6645
    gran.poly <- polygonGate(filterId = "Granulocytes",.gate=granulocytes@filter)
    cd66.45.poly <- polygonGate(filterId = "CD45+CD66- (Non Granulocytes)",.gate  =CD66negCD45pos@filter)
    return(list(gran.poly=gran.poly,cd66.45.poly=cd66.45.poly))
  })
  
  # ------ aggiugiamo granulociti e CD66-CD45+ al gating tree
  lista_poly_gran<-lapply(1:length(gs),function(i){
    list_poly_gran_cd45[[i]][[1]]
  })
  lista_poly_cd45<-lapply(1:length(gs),function(i){
    list_poly_gran_cd45[[i]][[2]]
  })
  names(lista_poly_gran)<-sampleNames(gs)
  # e aggiungiamolo al gating tree a partire dalle live cells.
  nodeID0<-add(gs,lista_poly_gran,parent="Live cells")
  # Creiamo il gateObject della pop Cd66-CD45+
  names(lista_poly_cd45)<-sampleNames(gs)
  # aggiungiamolo a partire sempre dalle Live cells
  nodeID0<-add(gs,lista_poly_cd45,parent="Live cells")
  recompute(gs) # modifichiamo effettivamente il gatingSet in base al nuovo albero
  # Estraiamo il flowSet dei granulociti.
  fs.gran <- getData(gs, "Granulocytes")
  fs.66n45p<-getData(gs,"CD45+CD66- (Non Granulocytes)")
  return(list(fs.gran=fs.gran,gs=gs,fs.66n45p=fs.66n45p,lista_poly_gran=lista_poly_gran,lista_poly_cd45=lista_poly_cd45))
}

# -----------------------------  Gating CD45+CD66- to obtain monocytes(HLADR+ CD14+) e HLADR+ CD14-  ---------------------------------------------------------
gating_ngran_to_DRpos<-function(i){
  # adesso partiamo dalla popolazione dei non Granulociti,ovvero CD45+Cd66-,prendiamo il loro flowSet dal GatingSet
  # ------- calcoliamo la HLADR-CD14- pop
  # Osserviamo questa popolazione rispetto al marker CD14: show(autoplot(fs.66n45p[[2]],"V500-A"))
  # come vediamo dal grafico ci sono due picchi vicini di cui uno piu' alto dell'altro.
  # bimodal = Logical. If TRUE, it returns a cutoff that splits population closer to 50-50, when
  # there are more than two peaks.
  # twin.factor = a value in [0,1] that is used to exclude twinpeaks
  # Calcoliamo il treshold piu' a destra della distribuzione(upper=T) rimuovendo eventuali picci gemelli (cioe' picchi piu' o meno allo stesso livello 
  # vicini,infatti twin.factor = .8) e cercando di dividere la popolazione 50/50(bimodal=t). Grazie a bimodal=T e twin.factor = 0.8 il threshold e' posizionatio fra
  # i due picchi invece che dopo i due picchi nonostante upper=T.
  lista_poly_mono_drp<-lapply(1:length(gs),function(i){
    #CD14.gate <- deGate(fs.66n45p[[i]], channel = c(fluorochrome.chans['CD14']), upper=T, bimodal = T, twin.factor = 0.8,alpha=.2,magnitude = .7)
    CD14.gate <- deGate(fs.66n45p[[i]], channel = c(fluorochrome.chans['CD14']), upper=T, alpha=0.9)
    
    # prendiamo la popolazione a sinistra di questo treshold quindi quella che comprende solo il picco piu' grande,escluendo il picco piu' piccolo.
    CD14neg <- flowDensity(fs.66n45p[[i]], channels = c(fluorochrome.chans['HLA-DR'], fluorochrome.chans['CD14']), position = c(NA, F), gates = c(NA, CD14.gate))
    # adesso da questa popolazione CD14- calcoliamo il treshold per HLA-DR marker(channel: BV605-A)
    # visualizziamo la densita' di questa pop CD14- rispetto al HLA-DR marker : show(autoplot(CD14neg@flow.frame,"BV605-A"))
    # c'e' un picco dominante piu' alto e altri picchi piu' piccoli e piatti vicino.
    # Calcoliamo tutti i tresholds possibili(all.cuts=T) 
    # Poi prendiamo un secondo threshold nell'estrema destra,dopo tutti i picchi.
    HLADR.gate <- c(deGate(CD14neg, channel = c(fluorochrome.chans['HLA-DR']), all.cuts = T, upper = T, percentile = NA),
                    deGate(CD14neg, channel = c(fluorochrome.chans['HLA-DR']), use.upper = T, upper = T))
    # calcoliamo la densita'  dei valori di espressione della popolazione Cd14- rispetto al marker  HLA-DR
    maxDens <- density(getflowFrame(CD14neg)@exprs[,c(fluorochrome.chans['HLA-DR'])])
    # prendiamo i valori di espressione da cui abbiamo calcolato la densita'(cioe' maxDens$x  ) sottriamo
    # 1.5,e troviamo l'indice del valore minimo,quindi il valore di espressione piu' basso della pop CD14- per HLA-DR
    start.idx <- which.min(abs(maxDens$x - 1.5))
    # con maxDens$y[start.idx:length(maxDens$y)] consideriamo tutti i valori di densita' di HLA-DR nella pop CD14-
    # prendiamo l' indice del valore massimo con  which.max(maxDens$y[start.idx:length(maxDens$y)])
    # e troviamo il valore di espressione che corrisponde alla densita' massima con maxDens$x[which.max(maxDens$y[start.idx:length(maxDens$y)])
    # Fra i threshold calcolati prima troviamo l' indice dei tresholds maggiori di questo valore di espressione sommatti al valore minimo di expr meno 1.
    # Di tutti treshold trovati prendiamo il treshold minimo.
    HLADR.gate <- min(HLADR.gate[which(HLADR.gate > maxDens$x[which.max(maxDens$y[start.idx:length(maxDens$y)]) + start.idx - 1])])
    # in pratica il nostro scopo e' quello di trovare un treshold piu' vicino possibile al picco piu' grande dalla parte destra. 
    # Quindi considerando tutti i tresholds piu' grandi del valore di espressione con densita' piu' grande,
    # stiamo di fatto prendendo tutti i thresholds a destra del picco grande di cui poi scegliamo quello 
    # piu' piccolo quindi quello piu' vicino al picco piu' grande. infatti per il primo flowFrame HLADR.gate vale 2.691278,quindi appena appena dopo il picco piu' grande.
    #  prendiamo la popolazione a sinistra di questo treshold,ovvero la pop HLA-DR-
    
    #HLADR.gate <- HLADR.gate*1.20
    
    DRn.14n <-flowDensity(CD14neg, channels = c(fluorochrome.chans['HLA-DR'], fluorochrome.chans['CD14']), position = c(F, NA), gates = c(HLADR.gate, NA))
    
    # ---- calcoliamo la monocites pop a partire da CD66-CD45+
    # notSubFrame = Remove a subset of a FlowFrame object specified by gates from the flowDensity method. It comes
    # in handy when one needs the complement of a cell population in the input flow cytometry data.
    # dunque rimuoviamo tutta questa popolazione (ovvero la pop CD14-HLADR-) dalla popolazione Cd45+CD66-,
    # osserviamo prima la pop CD66-CD45+ rispetto ai marker HLADR e CD14: show(autoplot(fs.66n45p[[1]],"BV605-A","V500-A"))
    temp <- notSubFrame(fs.66n45p[[i]], channels = c(fluorochrome.chans['HLA-DR'], fluorochrome.chans['CD14']), position = c(F, NA),filter= DRn.14n@filter)
    # mostriamo lo stesso grafico dopo la rimozione,ovvero plottiamo la pop rimanente: show(autoplot(temp@flow.frame,"BV605-A","V500-A"))
    # la pop rimanente,ovvere temp,e' la pop CD14+HLADR+
    # ruotiamo questa pop temp rimanente
    rot <- rotate.data(temp@flow.frame,c(fluorochrome.chans['HLA-DR'], fluorochrome.chans['CD14']),theta = -pi/4)
    # calcoliamo un treshold sulla pop ruotata.
    # visualizziamo questa pop temp ruotata rispetto a HLA-DR: show(autoplot(rot$data,"BV605-A"))
    
    # visualizziamo anche la pop temp ruotata rispetto a CD14: show(autoplot(rot$data,"V500-A"))
    #  due picchi di cui uno molto piu' grande dell'altro nel caso di HLA-DR
    # un unico picco grande nel caso di CD14
    # after.peack = Logical.  If TRUE, it returns a cutoff that is after the maximum peaks,  when there are more than two peaks
    # Questo threshold sta dopo il picco massimo, e in questo caso lo poniamo fra i due picchi. nel primo flowFrame e' pari a 0.1823185
    hladr.gate <- deGate(rot$data,fluorochrome.chans['HLA-DR'], upper=T,after.peak = T,twin.factor = .8,tinypeak.removal = 1/30)
    #Attenzione: da rimuovere nei PNG files.
    
    hladr.gate <- hladr.gate * 0.25
    
    # HLADRpos.CD14neg.flowD <-flowDensity(CD14neg.flowD, channels = c(fluorochrome.chans['HLA-DR'], fluorochrome.chans['CD14']), position = c(T, NA), gates = c(HLADR.gate, NA))
    # prendiamo la popolazione a sinistra di questo treshold (hladr.gate) per quanto riguarda HLA-DR, per quanto riguarda la densita' di CD14 per questa pop temp ruotata,
    # calcoliamo il threshold all'estrema destra (upper=c(NA,T),use.upper=c(F,T)) e prendiamo la popolazione a sinistra di questo treshold ( quindi in pratica prendiamo solo il picco grande, vedi il grafico di CD14)
    # ruotiamo la popolazione risultante della stessa quantita' di prima per annullare la rotazione precedente e otteniamo la popolazione dei monociti(HLADR+CD14+)
    monocytes <-  rotate.fd(flowDensity(rot$data, channels = c(fluorochrome.chans['HLA-DR'], fluorochrome.chans['CD14']), position = c(F, F), gates = c(hladr.gate, NA),upper=c(NA,T),use.upper=c(F,T)),angle = rot$theta)
    # Guardiamo la pop madre ovvero fs.66n45p senza monocytes,per capire dove originariamente erano piazzati i monocytes.
    # temp_prova <- notSubFrame(fs.66n45p[[1]], channels = c(fluorochrome.chans['HLA-DR'], fluorochrome.chans['CD14']), position = c(F, NA),filter= monocytes@filter)
    # show(autoplot(temp_prova@flow.frame,"BV605-A","V500-A"))
    # visualizziamo la pop originaria totale iniziale con il gate di questi monocytes
    #show(ggcyto(fs.66n45p[[1]],aes(x="BV605-A",y="V500-A")) + geom_hex() + geom_path(data=as.data.frame(monocytes@filter))) # sovrappongo due plots,in geom_path infatti gli dico di considerare un'altro dataset
    #indx <- which(monocytes@filter[,1]>=3.0)
    #monocytes@filter[indx,1]<- monocytes@filter[indx,1]*0.78
    #show(ggcyto(fs.66n45p[[1]],aes(x="BV605-A",y="V500-A")) + geom_hex() + geom_path(data=as.data.frame(monocytes@filter)))
    
    # ----- calcoliamo la pop HLADR+CD14-
    # prendiamo la popolazione a destra di questo threshold (hladr.gate) quindi la pop HLADR+, 
    # e ruotiamola della stessa quantita' di prima per annullare la rotazione precedente e otteniamo 
    # la pop HLADR+CD14-
    DRp.14n <-  rotate.fd(flowDensity(rot$data, channels = c(fluorochrome.chans['HLA-DR'], fluorochrome.chans['CD14']), position = c(T, NA), gates = c(hladr.gate, NA)),angle = rot$theta)
    # DRp.14n(output di rotate.fd) e' ancora un Cell population object con la matrice di espressione ruotata
    # visualizziamo la pop madre fs.66n45p  con il gate di questi DRp.14n
    # show(ggcyto(fs.66n45p[[1]],aes(x="BV605-A",y="V500-A")) + geom_hex() + geom_path(data=as.data.frame(DRp.14n@filter))) # sovrappongo due plots,in geom_path infatti gli dico di considerare un'altro dataset
    
    # ------- Descrizione delle 3 pop:DRn.14n,monocytes,DRp.14n
    # Dunque puoi vedere la differenza fra monocytes e DRp.14n:
    # - i monocytes sono cellule con valori di CD14(channel V500-A ) intorno al 3,
    # quindi molto piu' alti ma non vanno oltre il 3 per quanto riguarda l'HLADR(channel BV605-A).
    #  Dunque possiamo dire che i monocytes sono CD14+HLADR-. anche Se cmq alcuni hanno un valore alto di HLADR,quindi potremmo definirli anche CD14+HLADR+ 
    #- DRp.14n: hanno valori di HLADR intorno al 3,mentre CD14 si mantiene sotto il 3. Dunque per questo possiamo definirle CD14-HLADR+  
    #- DRn.14n hanno cellule sempre sotto il 2 sia per CD14 sia per HLADR e quindi sono  CD14-HLADR-
    #  in genere la pop neg sono quelle che rientrano nella prima decade(quindi minore di 1) ma non e' una regola.
    
    # dunque adesso abbiamo i boudaries della pop DRp.14n(appena trovata)
    # i boundaries della pop DRn.14n(trovati prima che avevamo rimosso allo scopo di ottenere meglio DRp.14n)
    # i boudaries della pop monocytes
    # questi boundaries potrebbero intersecarsi. Cerchiamo dunque di ottenere dei boundaries che non si intersecano.
    
    
    #------ confronto boundaries delle tre pops
    # creiamo i polygon objects dei boundaries di ognuna delle due pop.
    poly1 <- SpatialPolygons(list(Polygons(list(Polygon(DRp.14n@filter)), ID=c("c"))))                       
    poly2 <- SpatialPolygons(list(Polygons(list(Polygon(DRn.14n@filter)), ID=c("c"))))
    #gintersect = Function for testing if the geometries have at least one point in common or no points in common
    #if the two poligons have at least one point in common (quindi l'ouput di gIntersect e' True)
    if(gIntersects(poly1,poly2)){
      # gDifference = Function for determining the difference between the two given geometries.
      # Returns the regions of spgeom1 that are not within spgeom2
      n<-colnames(DRp.14n@filter) # nota che i nomi delle colonne di DRp.14n e DRn.14n sono le stesse
      diff.1<-rgeos::gDifference(poly2,poly1)
      # quindi diff.1 e' la regione di poly2 ovvero di DRn.14n non contenuta dentro DRp.14n.
      # riassegnamo le coordinate dei boundaries di DRn.14n come le coordinate della sua regione non contenuta in DRp.14n.
      DRn.14n@filter<- diff.1@polygons[[1]]@Polygons[[1]]@coords
      # riattacchiamo i nomi delle colonne.
      colnames( DRn.14n@filter)<-n
    }
    # Confrontiamo adesso la popolazione DRn.14n con la pop dei monocytes
    poly1 <- SpatialPolygons(list(Polygons(list(Polygon(monocytes@filter)), ID=c("c"))))
    poly2 <- SpatialPolygons(list(Polygons(list(Polygon(DRn.14n@filter)), ID=c("c"))))
    # facciamo praticamente la stessa cosa di prima.
    if (gIntersects(poly1,poly2))
    {
      n<-colnames(DRn.14n@filter)
      diff.1<-rgeos::gDifference(poly1,poly2)
      # modifichiamo dunque il filter originario dei monocytes con le coordinate del loro gate non contenuto
      # dentro il gate di  DRn.14n
      monocytes@filter<- diff.1@polygons[[1]]@Polygons[[1]]@coords
      colnames(monocytes@filter)<-n
    }
    # Ricalcoliamo i tre population objects a partire dalla popolazione madre indicando come boundaries dei loro gates i nuovi boundaries modificati.
    
    monocytes<-  flowDensity(fs.66n45p[[i]], channels = c(fluorochrome.chans['HLA-DR'], fluorochrome.chans['CD14']), position = c(T, F),filter=monocytes@filter)
    DRp.14n <-  flowDensity(fs.66n45p[[i]], channels = c(fluorochrome.chans['HLA-DR'], fluorochrome.chans['CD14']), position = c(T, F),filter=DRp.14n@filter)
    DRn.14n <-  flowDensity(fs.66n45p[[i]], channels = c(fluorochrome.chans['HLA-DR'], fluorochrome.chans['CD14']), position = c(T, F),filter=DRn.14n@filter)
    
    # ----- creiamo i polygon gates dei filters
    mono.poly <- polygonGate(filterId = "HLADR+ CD14+",.gate=monocytes@filter)
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
  
  
  
  #------ aggiungiamo la pop HLADR+CD14+(monocytes) e HLADR+ CD14- al gating tree
  # stessa cosa solita,creiamo il polygon objectGate object di ogni pop e aggiorniamo il gating tree
  names(list_mono.poly)<-sampleNames(gs)
  nodeID0<-add(gs,list_mono.poly,parent="CD45+CD66- (Non Granulocytes)")
  names(list_DRp.14n.poly)<-sampleNames(gs)
  nodeID0<-add(gs,list_DRp.14n.poly,parent="CD45+CD66- (Non Granulocytes)")
  names(list_DRn.14n.poly)<-sampleNames(gs)
  nodeID0<-add(gs,list_DRn.14n.poly,parent="CD45+CD66- (Non Granulocytes)")
  names(list_monocytes)<-sampleNames(gs)
  recompute(gs)
  fs.mono<-getData(gs,"HLADR+ CD14+")
  fs.DRp.14n<-getData(gs,"HLADR+ CD14-")
  fs.DRn.14n<-getData(gs,"HLADR- CD14-")
  return(list(fs.mono=fs.mono,fs.DRp.14n=fs.DRp.14n,fs.DRn.14n=fs.DRn.14n,gs=gs,list_monocytes=list_monocytes,list_DRp_14n=list_DRp_14n,list_DRn_14n=list_DRn_14n))
}


#------  Gating per arrivare ai Classical monocytes(CD16-) dai monocytes(HLADR+ CD14+) --------------
gating_mono_to_classic_mono<-function(i){
  # There are at least three types of monocytes in human blood:
  # - The classical monocyte is characterized by high level expression of 
  # the CD14 cell surface receptor (CD14++ CD16− monocyte)
  # - The non-classical monocyte shows low level expression of CD14 and additional co-expression 
  # of the CD16 receptor (CD14+CD16++ monocyte).[4]
  # - The intermediate monocyte with high level expression of CD14 and low level expression of CD16 (CD14++CD16+ monocytes).
  # Dunque I classical monocytes si individuano sfruttando il marker CD16 la cui espressione differisce fra i classical e i non classical monoCyte.
  # Dunque dovremmo partire dalla popolazione dei monocytes ovvero,HLADR+ CD14+ per trovare il gate dei CD16.
  # Tuttavia Mehrnoush ha scoperto che: "I found that finding the CD16 gate using the CD45+CD66-/CD14-CD16- population was easier than using the CD14+ 
  # population, because the CD14+ cells smear(si spalmano,strisciano) from CD16- to CD16+"
  # insomma ha visto che e' piu' facile calcolare il gate della CD16 sulla popolazione CD45+CD66- 
  # e e poi ricalcolandolo sulla popolazione CD14-CD16- rispetto che direttamente sulla popolazione CD14+
  list_classicalmono_poly<-lapply(1:length(gs),function(i){
    CD14.gate <- deGate(fs.66n45p[[i]], channel = c(fluorochrome.chans['CD14']), upper=T, bimodal = T, twin.factor = 0.8,alpha=.2,magnitude = .7)
    # Dunque consideriamo la popolazione madre ovvero 66n45p rispetto alla densita' di CD16 (channel FITC-A): show(autoplot(fs.66n45p[[3]],"FITC-A"))
    # Come vedi c'e' un unico grande picco centrale a sinistra (fra 1 e 2) quindi la maggior parte delle cellule e' negativa per CD16.
    # calcoliamo il threshold di CD16 sulla popolazione madre nella parte della coda del grafico (upper=T),quindi dopo il picco in questo caso,nel punto in cui cambia la slope.
    CD16temp.gate <- deGate(fs.66n45p[[i]], channel = c(fluorochrome.chans['CD16']), upper = T)
    # dunque prendiamo la popolazione a sinistra del threshold rispetto a CD16 (quindi la pop CD16-) 
    # e la pop a sinistra del treshold di CD14 che abbiamo calcolato prima(dunque una pop CD14-)
    # Dunque alla fine ho il Cell Population object di una pop CD14-CD16-
    #  visualizziamo la popolazione madre ovvero 66n45p rispetto alla densita' di CD14: show(autoplot(fs.66n45p[[1]],"V500-A","FITC-A"))
    CD14n.16n <- flowDensity(fs.66n45p[[i]], channels = c(fluorochrome.chans['CD14'], fluorochrome.chans['CD16']), position = c(F,F), gates = c(CD14.gate, CD16temp.gate))
    # calcoliamo il treshold di nuovo di CD16 stavolta sulla pop bivariata CD16-CD14-(use.upper=T e upper=T quindi dopo il picco) 
    # e prendiamo il minimo fra questo treshold e il threhold di CD16 calcolato sulla popolazione madre(ovvero fs.66n45p)
    # Questo treshold minimo e' il treshold di CD16 per trovare i classical monocytes.
    # Visualizziamo la pop CD14-CD16- rispetto a CD16: 
    # CD16monocyte.gate  <- min(deGate( CD14n.16n, channel = c(fluorochrome.chans['CD16']), use.upper = T, upper = T, alpha=0.05),
    #                           CD16temp.gate)
    all.peaks<-getPeaks(monocytes[[i]],channel = c(fluorochrome.chans['CD16']),tinypeak.removal = 0.1)
    if(max(all.peaks$Peaks)>=3.0){
      CD16monocyte.gate  <- deGate(monocytes[[i]], channel = c(fluorochrome.chans['CD16']), upper = T, alpha=0.9,tinypeak.removal = 0.1)
    }else{
      CD16monocyte.gate  <- deGate(monocytes[[i]], channel = c(fluorochrome.chans['CD16']), upper = T, alpha=0.9,tinypeak.removal = 0.1)*1.10
    }
    # Attenzione= da rimuovere nei files  PNG . Ho cambiato min con max e alpha=0.0001 da 0.05 a 0.0001 rimuovendo tinypeak.removal=0.02. Questo allo scopo di spostare il treshold piu' a destra.
    monocytes<-monocytes[[i]]
    #CD16monocyte.gate <- deGate(monocytes.flowD, channel = c(fluorochrome.chans['CD16']), use.upper = T, upper = T, tinypeak.removal = 0.2)
    # Calcoliamo dunque la popolazione dei classical monocytes a partire dai Monocytes usando il treshold di CD16 trovato sulla popolazione CD14-CD16- o sulla quella madre(in base a quale sia il minimo)
    # A me interessano infatti i valori a sinistra(quindi CD16-) rispetto a questo treshold.
    # rappresentiamo i monocytes rispetto a CD16: show(autoplot(monocytes[[1]]@flow.frame,"FITC-A"))
    classicalmono<- flowDensity(monocytes, channels = c(fluorochrome.chans['CD14'], fluorochrome.chans['CD16']), position = c(NA, F), gates = c(NA, CD16monocyte.gate))
    # usiamo setdiff() per trovare gli indici dei non classical monocytes
    # The elements of setdiff(x,y) are those elements in x but not in y.
    nonclassicalmono.index <- setdiff(monocytes@index, classicalmono@index)
    # Come al solito costruiamo il PolygonGate Object del filter dei classical monocytes e aggiungiamo il nodo al gating tree.
    classicalmono.poly <- polygonGate(filterId = "Classical Monocytes",.gate  =classicalmono@filter)
    return(list(classicalmono.poly=classicalmono.poly,CD16monocyte.gate=CD16monocyte.gate))
  })
  list_classicalmono_only_poly <- lapply(1:length(gs),function(i){
    list_classicalmono_poly[[i]]$classicalmono.poly
  })
  list_CD16mon_gate <- lapply(1:length(gs),function(i){
    list_classicalmono_poly[[i]]$CD16monocyte.gate
  })
  
  names(list_classicalmono_only_poly)<-sampleNames(gs)
  nodeID0<-add(gs,list_classicalmono_only_poly,parent="HLADR+ CD14+")
  recompute(gs)
  fs.classmono<- getData(gs,"Classical Monocytes")
  return(list(fs.classmono=fs.classmono,gs=gs,list_CD16mon_gate=list_CD16mon_gate))
}  

# versione alternativa che non mi piace
# gating_mono_to_classic_mono <- function(i){
#   Mono <- fsApply(fs.mono, function(x){
#     x.rot<-x
#     gate<- deGate(x.rot,fluorochrome.chans["CD14"],all.cuts = T,percentile=NA,sd.threshold = T)[1] 
#     gate.prctl<- deGate(x.rot,fluorochrome.chans["CD14"],use.percentile = T,percentile=.1)
#     temp <- notSubFrame(x.rot,c(fluorochrome.chans["CD14"],fluorochrome.chans["CD16"]),position = c(F,NA),gates=c(gate,NA))
#     cd16.gate.hi <-deGate(getflowFrame(temp),fluorochrome.chans["CD16"],upper=T,percentile=NA,tinypeak.removal = .9,alpha=.3,magnitude = .8)
#     n.mono<- flowDensity(x.rot,c(fluorochrome.chans["CD14"],fluorochrome.chans["CD16"]),position = c(T,T),gates=c(gate.prctl,cd16.gate.hi))
#     mono<- flowDensity(x.rot,c(fluorochrome.chans["CD14"],fluorochrome.chans["CD16"]),position = c(T,F),gates=c(gate,cd16.gate.hi))
#     plotDens(x,c(fluorochrome.chans["CD14"],fluorochrome.chans["CD16"]))
#     lines(mono@filter,type="l",lwd=2)
#     lines(n.mono@filter,type="l",lwd=2)
#     return(list(mono=mono,nmono=n.mono))
#   })
# }

#---------------------------- Gating HLADR+CD14- to Bcells,pDC,mDC -----------------------------------------------------
gating_drp_14n_to_B_pdc_mdc<-function(){
  list_poly_B_pdc_mdc<-lapply(1:length(gs),function(i){
    # adesso facciamo il gating dalle HLADR+CD14- per ottenere le Bcells,pDC e mDC populations.
    # prendiamo il flowSet del GatingSet del nodo HLADR+ CD14-
    # osserviamo la densità del marker CD11c di questa pop: show(autoplot(fs.drp.14n[[11]],"APC-A"))
    # show(autoplot(fs.drp.14n[[11]],"PE-Cy7-A"))
    # Due picchi di cui uno più grosso a sinistra e uno più piccolo a destra. 
    # Calcoliamo un treshold(opzioni di default) per CD11 della pop HLADR+ CD14-.
    CD11c.gate <- deGate(fs.drp.14n[[i]], channel = c(fluorochrome.chans['CD11c']),upper=T,alpha=0.01)*1.20
    # prendiamo una pop temporanea dalla pop HLADR+ CD14(markers CD123 e CD11c) a sinistra di questo treshold,quindi una pop CD11c-
    # if(CD11c.gate<1.5){
    #   CD11c.gate <- CD11c.gate*1.35
    # }
    #print(CD11c.gate)
    temp <- flowDensity(fs.drp.14n[[i]], channels = c(fluorochrome.chans['CD123'], fluorochrome.chans['CD11c']), position = c(NA, F), gates = c(NA, CD11c.gate))
    # rappresentiamo temp rispetto al marker CD123: show(autoplot(temp,"PE-Cy7-A"))
    # grande picco centrale, con qualche piccolo mini-picco qua e la'.
    #Da questa popolazione temporanea(quindi solo cellule CD11c-) calcoliamo il treshold per CD123 posizionandolo da qualche parte dopo il picco,
    # esattamente nell'ultimo picco della distribuzione senza considera i picchi piccoli se no mi va troppo a destra.
    #in generale non considero i picchi più piccoli(tinypeak.removal=0.2) per evitare che mi mette il treshold troppo lontano da dove voglio
    # al risultato aggiungiamo 0.2 per spostare il treshold un po' più a destra.
    # Ricorda:
    #  use.upper = T Logical.  If TRUE, forces to return the inflection point based on the first (last)
    # peak if upper=F (upper=T). default is F
    # upper = if TRUE, finds the change in the slope at the tail of the density curve, if FALSE,
    # finds it at the head.
    # Nota che la coda e la testa dipendono da dove si trova la media della distribuzione(valore più probabile della variabile),cioè se abbiamo un picco a destra ,la coda e' a sinistra e viceversa.
    # Quindi quando use.upper=T o =F estrare un preciso treshold nella coda(T) o nella testa(F)
    # Se use.upper=T e upper=F allora considera il primo picco della distribuzione
    # Se use.upper=T e upper=T considera l'ultimo picco della distribuzione 
    # Se use.upper=F torna il significato originale di upper
    # in questo caso ultimo picco senza considerare i  successivi picchi piccoli pari a 0.2 di proporzione (così il treshold sta un po' più a sinistra)
    # aggiungo inoltre 0.4,al posto di 0.2. Cancellare nei PNG files
    # show(autoplot(fs.drp.14n[[2]],"PE-Cy7-A"))
    CD123.gate <- deGate(temp, channel = c(fluorochrome.chans['CD123']), use.upper = T, upper = T, tinypeak.removal = 0.2,alpha = 0.001)*1.10
    print(identifier(fs.drp.14n[[i]]))
   
    # Calcoliamo inoltre il threshold posizionandolo da qualche parte a sinistra del picco,esattamente 
    # nel primo picco della distribuzione senza considerare i picchi piccoli per evitare che 
    # il treshold vada troppo a sinistra, inoltre sottraggo 0.05 per spostarlo solo un po' più a sinistra 
    #CD123.lo.gate <- deGate(temp, channel = c(fluorochrome.chans['CD123']), use.upper = T, upper = F, tinypeak.removal = 0.2) - 0.05 # primo picco della distribuzione senza considerare i precedenti picci piccoli(così il threshold parte da un po' più a destra)
    # prendiamo la popolazione mDC (CD123++,CD11c+) a partire dalla popolazione madre HLADR+CD14-, considerando le cellule a destra del treshold di CD123 sinistro(picco grande + picchi piccoli destri) 
    # e a destra del threshold Cd11c (picco piccolo a destra)
    mDC <- flowDensity(fs.drp.14n[[i]], channels = c(fluorochrome.chans['CD123'], fluorochrome.chans['CD11c']), position = c(F, T), gates = c(CD123.gate, CD11c.gate)) 
    # Attenzione da rimuovere nelle PNG files: al posto di position = c(F, T), gates = c(CD123.gate, CD11c.gate)) con  position = c(T, T), gates = c(CD123.low.gate, CD11c.gate)
    # prendiamo la popolazone pDC(CD123++) considerando le cellule con densità a destra del treshold CD123 destro, 
    #se vedi nel grafico corrispondente è la parte dopo l'unico picco grande,quindi comprende i i mini-picchi intorno al 3.
    # le pDC infatti sono cellule che esprimono molto CD123(i due + significano alta expressione) ma sono indifferenti nei confronti di CD11c.
    pDC<- flowDensity(temp, channels = c(fluorochrome.chans['CD123'], fluorochrome.chans['CD11c']), position = c(T, NA), gates = c(CD123.gate, NA))
    # prendiamo la popolazione Bcells (CD123+) considerando le cellule a sinistra del treshold  CD123 destro,quindi prendiamo il picco grande insieme ai picchi piccoli sinistri
    Bcells<- flowDensity(temp, channels = c(fluorochrome.chans['CD123'], fluorochrome.chans['CD11c']), position = c(F, NA), gates = c(CD123.gate, NA))
    # da popolazione poi prendiamo le cellule a destra del treshold sinistro CD123 in modo da prende solo le cellule con densità in corrispondenza del picco centrale.
    # ATTENZIONE: nei PNG files è position T con cd123.low.gate nella seconda Bcells che qui è cancellata
    #Bcells<- flowDensity(Bcells, channels = c(fluorochrome.chans['CD123'], fluorochrome.chans['CD11c']), position = c(F, NA), gates = c(CD123.gate, NA)) 
    # Dunque le B cells sono cellule che esprimono CD123 ma non troppo,cioè non a livelli comparibili con le altre due pop.
    mdc.poly <- polygonGate(filterId = "mDC",.gate  =mDC@filter)
    pdc.poly <- polygonGate(filterId = "pDC",.gate  =pDC@filter)
    bcell.poly <- polygonGate(filterId = "B cells",.gate  =Bcells@filter)
    return(list(mdc.poly=mdc.poly,pdc.poly=pdc.poly,bcell.poly=bcell.poly))
  })
  lista_mdc_poly <- lapply(1:length(list_poly_B_pdc_mdc),function(i){
    list_poly_B_pdc_mdc[[i]][[1]]
  })
  lista_pdc_poly <- lapply(1:length(list_poly_B_pdc_mdc),function(i){
    list_poly_B_pdc_mdc[[i]][[2]]
  })
  lista_bcell_poly <- lapply(1:length(list_poly_B_pdc_mdc),function(i){
    list_poly_B_pdc_mdc[[i]][[3]]
  })
  
  #------ aggiunta gates al gating tree
  #A questo punto facciamo le solite operazioni. Creiamo il polygonGate del filter object di ogni popolazione
  # e poi aggiungiamo il polygon gate al gating tree. la popolazione madre di tutte e tre le popolazione e' sempre
  # la popo HLADR+ CD14-.
  names(lista_mdc_poly)<-sampleNames(gs)
  nodeID0<-add(gs,lista_mdc_poly,parent="HLADR+ CD14-")
  names(lista_pdc_poly)<-sampleNames(gs)
  nodeID0<-add(gs,lista_pdc_poly,parent="HLADR+ CD14-")
  names(lista_bcell_poly)<-sampleNames(gs)
  nodeID0<-add(gs,lista_bcell_poly,parent="HLADR+ CD14-")
  recompute(gs)
  return(list(gs=gs,lista_mdc_poly=lista_mdc_poly,lista_pdc_poly=lista_pdc_poly,lista_bcell_poly=lista_bcell_poly))
}


#-------------- Gating HLADR-CD14- to get T cells(CD3+), gdT(CD3+gd+) e CD3-  --------------------------------------------------------
gating_drn_14n_to_T_gdt_cd3<-function(){
  list_t_gdt_cd3_poly<-lapply(1:length(gs), function(i){
    # adesso ripartiamo dalla pop HLADR-CD14- ,ovvero DRn.14n, figlia di CD45+CD66-, che non abbiamo ancora aggiunto al gating tree.
    # Da questa pop  otteniamo la pop gdT.
    # Visualizziamo questa pop rispetto alla densità di CD3: show(autoplot(fs.drn.14n[[1]],"PE-CF594-A"))
    # Ci due picchi uno molto alto alla'estrema destra e uno basso all'estrema sinistra.
    # Calcoliamo il treshold di default di CD3. Che si trova in mezzo ai due picchi.
    CD3.gate <- deGate(fs.drn.14n[[i]], channel = c(fluorochrome.chans['CD3']),upper=F,alpha = 0.05)*0.90 
    #print(CD3.gate)
    #  prendiamo la pop CD3Tcells(CD3+) considerando le cellule con densità a destra di questo treshold. 
    CD3Tcells<- flowDensity(fs.drn.14n[[i]], channels = c(fluorochrome.chans['CD3'], scat.chans['SSC-A']), position = c(T, NA), gates = c(CD3.gate, NA))
    # prendiamo la pop CD3neg (CD3-) considerando le cellule con densità a sinistra di questo treshold.
    CD3neg <- flowDensity(fs.drn.14n[[i]], channels = c(fluorochrome.chans['CD3'], scat.chans['SSC-A']), position = c(F, NA), gates = c(CD3.gate, NA))
    # visualizziamo la pop CD3Tcells(CD3+) rispetto al marker gd:show(autoplot(CD3Tcells@flow.frame,"PE-A"))
    # unico picco grande fra 1 e 2,qualche mini-picco ai lati.
    
    # Ricorda:
    #  use.upper = T Logical.  If TRUE, forces to return the inflection point based on the first (last)
    # peak if upper=F (upper=T). default is F
    # upper = if TRUE, finds the change in the slope at the tail of the density curve, if FALSE,
    # finds it at the head.
    # Nota che la coda e la testa dipendono da dove si trova la media della distribuzione(valore più probabile della variabile),cioè se abbiamo un picco a destra ,la coda e' a sinistra e viceversa.
    # Quindi quando use.upper=T o =F estrare un preciso treshold nella coda(T) o nella testa(F)
    # Se use.upper=T e upper=F allora considera il primo picco della distribuzione
    # Se use.upper=T e upper=T considera l'ultimo picco della distribuzione 
    # Se use.upper=F torna il significato originale di upper
    
    # calcoliamo il threshold sulla pop CD3T per il channel gd. Calcoliamo prima tutti i tresholds(all.cuts = T) nella coda della distribuzione(in questo caso si trova prima del picco grande)
    # rimuovendo quelli davvero troppo piccoli(0.001)
    # quindi calcoliamo il minimo di questi tresholds perchè vogliamo il treshold il piu' possibile spostato a sinistra.
    # gd.gate <- min(deGate(CD3Tcells, channel = c(fluorochrome.chans['gd']), all.cuts = T, upper = T, percentile = NA, tinypeak.removal = 0.001, alpha = 0.05))
    # # valore di gd.gate alla fine è 0.4086271 per il primo flowFrame.
    # # Consideriamo temporaneamente la popolazione CD3Tcells(CD3+) rispetto ai marker CD3 e gd.
    # temp <- getflowFrame(CD3Tcells)@exprs[, c(fluorochrome.chans['CD3'], fluorochrome.chans['gd'])]
    # ricalcoliamo il treshold per gd posizionandolo all'ultimo picco (esclusi picchi piccoli).
    gd.gate <-deGate(CD3Tcells@flow.frame,fluorochrome.chans['gd'],use.upper=T,upper=T,tinypeak.removal=.99,alpha=.05)
    # prendiamo una popolazione temporanea  gd+ dalla madre CD3T considerando le cellule con densità a destra di quest'ultimo treshold.
    temp1<- flowDensity(CD3Tcells@flow.frame, c(fluorochrome.chans['CD3'], fluorochrome.chans['gd']),position = c(NA,T),gates=c(NA,gd.gate))
    # Dunque temp1 e' già una popolazione CD3+ (perchè deriva da CD3T) ed è anche gd+ come abbiamo appena visto.
    # però io voglio le cellule di temp1 ancora più positive,con ancora più valori di espressione di CD3.
    # Dunque visualizziamo temp1 rispetto a CD3: show(autoplot(temp1@flow.frame,"PE-CF594-A"))
    # grande picco a destra,serie di piccoli picchi a sinistra
    # Su questa pop temp1 calcoliamo un threshold posizionandolo come primo picco della distribuzione(quindi a sinistra) senza considerare i picchi piccoli.
    cd3.gate.lo <- deGate (getflowFrame(temp1),fluorochrome.chans['CD3'],upper=F,use.upper=T,tinypeak.removal = .9,alpha=.1)
    # in questo modo ho ottenuto un treshold che mi consente di ottenere le cellule ancora più CD3+ del CD3.gate iniziale ( ma tenendo in considerazione che devono essere cmq cell gd+ infatti calcolo questo treshold su temp1)
    # Dunque riprendiamo la popolazione CD3Tcells(che sono solo cellule CD3+)
    # da questa pop prendiamo la popolazione delle gdTcells(CD3+gd+) considerando le cellule a destra del thresold di gd e del secondo treshold di CD3 trovato poco che mi consente di avere cellule ancora più positive per CD3.
    gdTcells<-flowDensity(CD3Tcells@flow.frame, c(fluorochrome.chans['CD3'], fluorochrome.chans['gd']),position = c(T,T),gates=c(cd3.gate.lo,gd.gate))
    
    #-----  creiamo i polygon gate objects
    tcell.poly <- polygonGate(filterId = "CD3+ T cells",.gate  = CD3Tcells@filter)
    cd3n.poly <- polygonGate(filterId = "CD3-",.gate  =CD3neg@filter)
    gd.poly <- polygonGate(filterId = "gd T cells",.gate  =gdTcells@filter)
    return(list(tcell.poly=tcell.poly,cd3n.poly=cd3n.poly,gd.poly=gd.poly,CD3.gate=CD3.gate))
  })
  list_t_poly<-lapply(1:length(gs),function(i){
    list_t_gdt_cd3_poly[[i]]$tcell.poly
  })
  list_gdt_poly<-lapply(1:length(gs),function(i){
    list_t_gdt_cd3_poly[[i]]$gd.poly
  })
  list_cd3_poly<-lapply(1:length(gs),function(i){
    list_t_gdt_cd3_poly[[i]]$cd3n.poly
  })
  list_CD3_gate<-lapply(1:length(gs),function(i){
    list_t_gdt_cd3_poly[[i]]$CD3.gate
  })
  
  #----- add gates to the gating tree
  # Facciamo le solite operazioni cioè creiamo il poligonGate object dei boundaries della popolazione e
  # aggiungiamo i gates al gating tree.
  # la popolazione T di cellule sono semplicemente cellule CD3+ indifferenti rispetto a gd.
  names(list_t_poly)<-sampleNames(gs)
  nodeID0<-add(gs,list_t_poly,parent="HLADR- CD14-")
  # la popolazione Cd3neg è la pop CD3- come dice lo stesso nome,indifferenti rispetto a gd.
  names(list_gdt_poly)<-sampleNames(gs)
  nodeID0<-add(gs,list_gdt_poly,parent="HLADR- CD14-")
  recompute(gs)
  # la popolazione gdT sono le T cells che sono anche positive al marker gd+ (quindi sono CD3+gd+ cells)
  names(list_cd3_poly)<-sampleNames(gs)
  nodeID0<-add(gs,list_cd3_poly,parent="HLADR- CD14-")
  recompute(gs)
  fs.tcell<- getData(gs,"CD3+ T cells")
  fs.cd3n<-getData(gs,"CD3-")
  return(list(fs.tcell=fs.tcell,fs.cd3n=fs.cd3n,gs=gs,list_gdt_poly=list_gdt_poly,list_CD3_gate=list_CD3_gate))
  
}


# cd45.gate <-fsApply(fs.live,deGate,channels.ind[4],upper=F,use.upper=T,percentile=NA,tinypeak.removal=.5,alpha=.01)*.95
# names(cd45.gate)<-sampleNames(fs.live)
# Lymph <- fsApply(fs.live, function(x){
#   
#   cd45.pos <- flowDensity(x,c(4,channels.ind[4]),position = c(NA,T),gates=c(NA,cd45.gate[identifier(x)]))
#   ss.gate <- deGate(getflowFrame(cd45.pos),4)
#   lymph <- flowDensity(x,c(4,channels.ind[4]),position = c(F,T),gates=c(ss.gate,cd45.gate[identifier(x)]))
#   gran <- flowDensity(x,c(4,channels.ind[4]),position = c(T,T),gates=c(ss.gate,cd45.gate[identifier(x)]))
#   plotDens(x,c(4,channels.ind[4]))
#   lines(lymph@filter,type="l",lwd=2)
#   lines(cd45.pos@filter,type="l",lwd=2)
#   return(list(lymph=lymph, gran=gran,cd45=cd45.pos))
# })
# fs.lymph <- lapply(Lymph, function(x) getflowFrame(x$lymph))
# fs.lymph  <- as(object=fs.lymph , Class="flowSet")


# -------------------------- Gating NKT and NK -------------------------------------------------------------


# Dunque un po' di teoria:
# Le NK cells sono cellule che esprimono alti livelli di CD16 e CD56,quindi sono cellule CD16+CD56+.
# a sua volta le NK possono differenziarsi in base al livello di CD56 come CD56dim e CD56 bright.
# C'è una pop di cellule che esprime alte concentrazioni di CD56(come le NK) ma anche di CD3 (come le T cells)
# quidni sono chiamate NKT cells,che sono CD56+CD3+
# Dunque abbiamo:
# CD56 + CD3 -=  NK cells, 
# CD56 - CD3 + = T lymphocytes, che abbiamo trovato nel precedente gating
# CD3 + CD56 + = NKT 
# A sua volta NK cells si dividono in molte subpopolazioni in base all'espressione di CD16,CD56,CD8,CD62L,CD38
# Inoltre vedi schema flowCytometry_theory_and_data_formats per conoscere cosa vogliono dire i simboli accanto ai markers.
# le NK cells ( e quindi anche le NKT) presentano inoltre l'espressione di CD16(quindi sono CD16+)
# Tuttavia qualche volta le NK presentano cellule CD16-,il motivo sarà sicuramente biologico ma è sconosciuto(osservato da rym e merhnoush su almeno due sample)



gating_nkt_nk<-function(){
  list_all_nkt_related_pops<-lapply(1:length(gs),function(i){
    # ------------ Per il momento cerchiamo di ottenere le NKT cells 
    # quindi partiamo dalla madre CD3+.
    # Dalle CD3+ cells(CD3T cells) facciamo il gating per ottenere le NKT cells,una subpopulazione delle T cells.
    # prendiamo il flowSet del GatingSet della CD3T population
    # visualizziamo questa pop in base alla densità di CD16: show(autoplot(fs.tcell[[6]],"FITC-A","BV650-A"))
    # grande picco a sinistra,un po' dopo l'1,mini picchi a destra di esso.
    # calcoliamo il treshold sulla pop madre CD3T cells per CD16 come ultimo picco della della distribuzione(semza considera i picchi piccoli successivi)
    # 1.8 circa quindi a destra del picco centrale,che e' un picco di negativi,perche' è intorno a 1.
    CD16NKT.gate <- deGate(fs.tcell[[i]], channel = c(fluorochrome.chans['CD16']), use.upper = T, upper = T, alpha = 0.05, tinypeak.removal = 0.4)
    # Ricorda:
    #  use.upper = T Logical.  If TRUE, forces to return the inflection point based on the first (last)
    # peak if upper=F (upper=T). default is F
    # upper = if TRUE, finds the change in the slope at the tail of the density curve, if FALSE,
    # finds it at the head.
    # Nota che la coda e la testa dipendono da dove si trova la media della distribuzione(valore più probabile della variabile),cioè se abbiamo un picco a destra ,la coda e' a sinistra e viceversa.
    # Quindi quando use.upper=T o =F estrare un preciso treshold nella coda(T) o nella testa(F)
    # Se use.upper=T e upper=F allora considera il primo picco della distribuzione
    # Se use.upper=T e upper=T considera l'ultimo picco della distribuzione 
    # Se use.upper=F torna il significato originale di upper
    
    # Visualizziamo la pop delle CD3T rispetto a CD56: show(autoplot(fs.tcell[[1]],"BV650-A"))
    # grande picco centrale intorno al 2,mini-picci a sinistra e a destra,coda a sinistra.
    # Intorno a 2 non e' una grande valori di espressione,dunque i cd56 positivi sono quelli intorno al 3(che sono molto poche)
    # getPeaks() = Find all peaks in density along with their indices
    # Calcoliamo dunque tutti i picchi di densità di CD56 nella pop CD3T(senza considerare i picchi piccoli)
    peaks <-getPeaks(fs.tcell[[i]], channel = c(fluorochrome.chans['CD56']), tinypeak.removal = 0.05)
    # calcoliamo il treshold nella pop CD3T per CD56. Prima calcoliamo il treshold di default(senza considerare i picchi piccoli) che e' pari a 0.8172373
    # poi calcoliamo il treshold come ultimo picco della distribuzione,quindi dopo il picco (quindi pari a 2.506033)
    # e poi calcoliamo il minimo fra i due quindi il risultato e' 0.8172373 ovvero prima del picco grande.
    CD56NKT.gate <-  min(deGate(fs.tcell[[i]], channel = c(fluorochrome.chans['CD56']), tinypeak.removal = 0.05),
                         deGate(fs.tcell[[i]], channel = c(fluorochrome.chans['CD56']), use.upper = T, upper = T, tinypeak.removal = 0.9))
    # Se questo treshold è minore del valore del picco che ha la massima altezza (quindi del picco grande intorno a 2)
    # non va bene perchè io voglio il threshold che mi separi i positivi dai negativi se e' prima del picco significa che è in mezzo ai negativi(il picco grande intorno a 2 e' un picco di negativi)
    # aggiungo inoltre un altro ifstatement perche' in un sample ha dato errore dato che CD56NKT.gate<peaks$Peaks[which.max(peaks$P.h) risultava NA
    if(is.na(CD56NKT.gate<peaks$Peaks[which.max(peaks$P.h)])==FALSE){
      if(CD56NKT.gate<peaks$Peaks[which.max(peaks$P.h)]){
        # ricalcoliamo il treshold considerandolo come(inflection point) ultimo picco della distribuzione. Per il primo flowFrame e' pari a  2.33,quindi a destra del picco grande.
        CD56NKT.gate  <- deGate(fs.tcell[[i]], channel = c(fluorochrome.chans['CD56']), use.upper = T, upper = T, tinypeak.removal = 0.9,alpha=.5)
      } 
    }
    # prendiamo la popolazione di CD3T considerando le cellule con densità a destra del treshold di CD16 e a destra del treshold di CD56. Quindi la pop CD16+CD56+. 
    CD3.16p56p <- flowDensity(fs.tcell[[i]], channels = c(fluorochrome.chans['CD16'], fluorochrome.chans['CD56']), position = c(T, T), gates = c(CD16NKT.gate, CD56NKT.gate))
    # prendiamo la popolazione a sinistra del treshold di CD16 e a destra del  treshold di CD56,quindi la pop CD16-CD56+
    CD3.1np56p <- flowDensity(fs.tcell[[i]], channels = c(fluorochrome.chans['CD16'], fluorochrome.chans['CD56']), position = c(F, T), gates = c(CD16NKT.gate, CD56NKT.gate))
    # prendiamo la popolazione a sinistra di CD16 e a sinistra del treshold di CD56,quindi la pop CD16-CD56-.
    CD3.16n56n <- flowDensity(fs.tcell[[i]], channels = c(fluorochrome.chans['CD16'], fluorochrome.chans['CD56']), position = c(F, F), gates = c(CD16NKT.gate, CD56NKT.gate))
    return(list(CD3.16p56p=CD3.16p56p,CD3.1np56p=CD3.1np56p,CD3.16n56n=CD3.16n56n,CD56NKT.gate=CD56NKT.gate,CD16NKT.gate=CD16NKT.gate))
  })
  
  #------- Adesso cerchiamo le NK cells,dunque partiamo dalla madre CD3-
  # NK cell gating: Discussion from Rym (email 5/15/2017)
  # For the NK cells: It is not unusual to see the NK cells shifting on the CD16 axis, we encountered this before,
  # we think it's biological but the reasons are unknown. Two of the individuals run for HVP have their NK cells CD16-,
  # CD56 dim. In that case, and since it's ONE distinct population, we move the quadrant gate to the right on the Y
  # axis (to not have the population split in two).
  # NK cells: The way we usually look at those gates when comes the time for analysis, we take all the counts/ MFI or
  # % from the 3 seperate gates : CD56 high NK, CD56 Dim CD16 +/- NK and CD56 -CD16+ NK.
  
  # Nota che fra myeloid cells e PNG cells(a cui si riferisce la gating Strategy) c'e' una differenza nella densita' di NK,infatti se confronti i plotdens,noti che la popolazione
  # in alto a destra nelle meyloid cells e' assente invece nei PNG e' presente,quindi devi adattare il codice a questa situazione.
  list_all_nk_pop<-lapply(1:length(gs),function(i){
    # prendiamo il flowFrame della pop CD3-
    # f.3n <-getflowFrame(CD3neg)
    # mia versione:
    f.3n<-fs.cd3n[[i]]
    # visualizziamo questo flowFrame in base a CD6 e CD56: show(autoplot(fs.cd3n[[14]],"FITC-A","BV650-A")),show(autoplot(fs.cd3n[[2]],"FITC-A",))
    # rimuoviamo le marginal cells da questa pop CD3 eliminando anche le negative cells con negatività inferiore a 0.1,
    x<-removeMargins(f.3n,fluorochrome.chans[c('CD16','CD56')],neg=.1);
    # x e' il flowFrame object con la matrice di espressione filtrata,ovvero senza più i marginal events.
    # visualizziamo questa pop filtrata in base a CD16 e CD56: show(autoplot(x,"FITC-A","BV650-A"))
    # prendiamo una pop temporanea da questa pop filtrata considerando le cellule a sinistra del treshold di CD56 moltiplicato per 0.95,
    # quindi spostiamo il threshold un po' più a sotto(in questo caso CD56 è nell'asse y) da 2.33 a 2.2135 per assicurarci di prendere solo il picco sotto, e quindi una pop CD3-C56-,indiferrente per CD16 per ora
    CD56NKT.gate <- list_all_nkt_related_pops[[i]]$CD56NKT.gate
    temp <- flowDensity(x,fluorochrome.chans[c('CD16','CD56')],position = c(NA,F),gates=c(NA,CD56NKT.gate*.95))
    # da questa pop CD3-CD56- calcoliamo il treshold per CD16
    # visualizziamo questa pop temp rispetto a CD16: show(autoplot(temp,"FITC-A"))
    # due picchi,uno grande a destra dopo l'1 e uno picco a sinistra vicino al 3.
    # calcoliamo il treshold per CD16 nella coda della distribuzione (il picco grande,ovvero la testa è a sinistra quindi la coda è a destra)
    # e infatti il threshold è dopo il picco (2.448732)
    cd16.gate.hi<-deGate(getflowFrame(temp),fluorochrome.chans['CD16'],upper=T,alpha=.8,tinypeak.removal = 0.5)*1.15
    # dalla pop x(CD3-) prendiamo la pop CD56+,quindi viene CD3-CD56+ considerando le cellule a destra (quindi sopra) il treshold per CD56
    temp <- flowDensity(x,fluorochrome.chans[c('CD16','CD56')],position = c(NA,T),gates=c(NA,CD56NKT.gate*.9))
    # un picco grande vicini al 2 con uno shoulder a sinistra.
    # calcoliamo il treshold  per CD16 da questa pop CD3-CD56+ senza considerare i picchi piccoli pari a 1/40 cercando di dividere la pop 50/50 (infatti bimodal = True) e
    # lo calcoliamo in base all deviazione standard con n.sd pari a 1.55
    # sd.threshold = if TRUE, uses ’n.sd’ times standard deviation as the threshold. Default value is set to ’FALSE’.
    # n.sd = an integer coefficient for the standard deviation to determine the threshold based on the standard deviation if ’sd.threshold’ is TRUE.
    # esso per il primo flowFrame e' pari a 2.15852 quindi dopo il picco grande.
    #cd16.gate.lo<- deGate(getflowFrame(temp),fluorochrome.chans['CD16'],tinypeak.removal = 1/40,sd.threshold=T,n.sd=1.55,bimodal=T,alpha=.5)
    print(identifier(f.3n))
    all.peaks <- getPeaks(getflowFrame(temp),fluorochrome.chans['CD16'],tinypeak.removal = 0.1)
    #codice da rimuovere nei PNG files
    if(length(all.peaks$Peaks)<2){
      if(all.peaks$Peaks<2.0){
        cd16.gate.lo <- deGate(getflowFrame(temp),fluorochrome.chans['CD16'],use.upper=T,upper=T,alpha=0.3,tinypeak.removal = 0.1)*0.95
      } else{
        cd16.gate.lo <- deGate(getflowFrame(temp),fluorochrome.chans['CD16'],upper=F,alpha=0.2,tinypeak.removal = 0.1)
      }
    }else if(length(all.peaks$Peaks)>=2){
      if(max(all.peaks$Peaks)<2.0){
        cd16.gate.lo <- deGate(getflowFrame(temp),fluorochrome.chans['CD16'],use.upper=T,upper=T,alpha=0.3)*0.90
      } else{
        cd16.gate.lo <- deGate(getflowFrame(temp),fluorochrome.chans['CD16'],upper=F,alpha=0.3,tinypeak.removal = 0.1)
      }
    }
    
    #print(all.peaks)
    # print(cd16.gate.lo)
    # dalla popolazione CD3-,prendiamo la popolazione a destra del treshold più alto di CD16  e a sinistra del treshold di CD56 (sotto perchè CD56 è nell'asse y)
    quad.2<-flowDensity(f.3n,fluorochrome.chans[c('CD16','CD56')],position = c(T,F),gates=c(cd16.gate.hi,CD56NKT.gate*.9))
    # la chiamiamo quad 2 perchè è la pop che corrisponde con il secondo quadrante,ovvero quell in basso a destra.
    # visualizziamo questa pop rispetto a CD56: show(autoplot(quad.2@flow.frame,"BV650-A"))
    # picco grande fra 1 e 1.5 ( ricorda il quad2 contiene le cells con valori di CD56 inferiori a 2.25 circa)
    # da questa pop CD16++CD56- calcoliamo il treshold per CD56 come ultimo picco(inflection point) della distribuzione
    # count.lim = minimum limit for events count in order to calculate the threshold.  Default is 20, returning NA as threshold.
    cd56.gate.lo<-deGate(getflowFrame(quad.2),fluorochrome.chans['CD56'],use.upper=T,upper=T,tinypeak.removal = .9,count.lim = 5)
    # esso è pari a 2.077995,quindi dopo il picco grande
    # se e' cd56.gate.lo e' inf allora settiamo cd56.gate come 10,se no prendiamo cd56.gate.lo così come.
    # a questo il tresholld di  cd56 finale è il minimo fra cd56.gate.lo e CD56NKT.gate*.93 (rimuoviamo gli NA)
    cd56 <- min(ifelse(is.infinite(cd56.gate.lo),yes = 10,no=cd56.gate.lo),CD56NKT.gate*.93,na.rm = F)
    # Dall popolazione madre CD3- prendiamo la popolazione a destra del CD16 più alto che abbiamo trovato e a sinistra di quest'ultimo treshold di CD56
    # Dunque otteniamo la pop CD56-CD16+ (che è anche CD3-)
    nk.56n16p <-flowDensity(f.3n,fluorochrome.chans[c('CD16','CD56')],position = c(T,F),gates=c(cd16.gate.hi,cd56))
    # prendiamo stavolta la pop a sinistra del CD16 più alto e dell'ultimo treshold di CD56,quindi
    # otteniamo la pop CD56-CD16-
    nk.56n16n<-flowDensity(f.3n,fluorochrome.chans[c('CD16','CD56')],position = c(F,F),gates=c(cd16.gate.hi,cd56))
    # prendiamo la pop temporenea NK con cellule con densità a sinistra del CD16 più basso e a destra dell'ultimo CD56 calcolato
    # quindi una pop CD16- e CD56+
    NK.temp<-flowDensity(f.3n,fluorochrome.chans[c('CD16','CD56')],position = c(F,T),gates=c(cd16.gate.lo,cd56))
    
    # infine prendiamo la popolazione a destra di CD16 più basso e a destra di ultimo CD56.
    # otteniamo la pop CD56+CD16+. Comprensde tutti i valori positivi di Cd56 a prescindere che siano dim,mod o high.
    nk.56p16p<-flowDensity(f.3n,fluorochrome.chans[c('CD16','CD56')],position = c(T,T),gates=c(cd16.gate.lo,cd56))
    # la differenza fra con la pop CD16+ di prima è che quest'ultimo comprende anche valori di fluorescenza moderate,mentre invece la prima CD16+ comprende solo bright values,
    # infatti è stata ottenuta considernado la destra del valore più alto del treshold i CD16
    # adesso visualizziamo la Nktemp rispetto a CD56: show(autoplot(NK.temp@flow.frame,"BV650-A"))
    # due picchi uno grande fra 2.5 e 3.0 e uno piccolo a 3.5
    # adesso calcoliamo il treshold sulla NK temp(CD16-Cd56+) in base alla  deviazione standard pari a 1 senza considerare i picchi piccoli
    cd56.gate.sd <- deGate(obj = getflowFrame( NK.temp),channel = fluorochrome.chans['CD56'],percentile=NA,sd.threshold=T,n.sd=1,tinypeak.removal = .8)
    # pari a 3.119
    cd56.gate.hi <- tail(getPeaks(getflowFrame( NK.temp),fluorochrome.chans['CD56'],tinypeak.removal=.1)$Peaks,1)-.8*sd(exprs(getflowFrame(NK.temp))[,fluorochrome.chans['CD56']])
    # print(identifier(f.3n))
    # print(cd56.gate.hi)
    # pari a 3.20
    cd56.gate.hi <- max(cd56.gate.sd,cd56.gate.hi)
    # print(cd56.gate.hi)
    # dunque 3.20,un treshold a destra del picco grande.
    # prendiamo la popolazione che ha CD16- e CD56+ ma solo leggermente positivo,positivi di poco diciamo. Dunque prendiamo la pop a sinistra(sotto)
    # il treshold cd56.gate.hi appena calcolato,ovvero circa 3
    nk.56dim16n<- flowDensity( NK.temp,fluorochrome.chans[c('CD16','CD56')],position = c(NA,F),gates=c(NA,cd56.gate.hi))
    # prendiamo la popolazione CD56+ ma stavolta con alta espressione di CD56(bright expression),quindi la pop a destra del treshold cd56.gate.hi,dopo il 3 quindi verso il picco di 3.5
    nk.56hi <-flowDensity( NK.temp,fluorochrome.chans[c('CD16','CD56')],position = c(NA,T),gates=c(NA,cd56.gate.hi))
    
    # costruiamo i polygon gates object
    nk.56hi.poly <- polygonGate(filterId = "CD56 Hi NK",.gate  =nk.56hi@filter)
    nk.56n16n.poly <- polygonGate(filterId = "CD56-CD16-",.gate  =nk.56n16n@filter)
    nk.56n16p.poly <- polygonGate(filterId = "CD56-CD16+ NK",.gate  =nk.56n16p@filter)
    nk.56dim16n.poly <- polygonGate(filterId = "CD56 dim CD16- NK",.gate  =nk.56dim16n@filter)
    nk.56p16p.poly <- polygonGate(filterId = "CD56 dim CD16+ NK",.gate  =nk.56p16p@filter)
    
    
    
    return(list(nk.56hi.poly=nk.56hi.poly,nk.56n16n.poly=nk.56n16n.poly,
                nk.56n16p.poly=nk.56n16p.poly,nk.56dim16n.poly=nk.56dim16n.poly,
                nk.56p16p.poly=nk.56p16p.poly,nk.56n16n=nk.56n16n,
                NK.temp=NK.temp,nk.56hi=nk.56hi))
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
  
  #--------- aggiugiamo gated pop NK(quindi dalla madre CD3-) al gating tree
  # stessa cosa che facciamo di solito creiamo i polygon gate dei filter e li aggiungiamo al gating tree.
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
  recompute(gs)
  fs.nk.56n16n<-getData(gs,"CD56-CD16-")
  return(list(fs.nk.56n16n=fs.nk.56n16n,gs=gs,list_all_nkt_related_pops=list_all_nkt_related_pops,
              list_nk_56n16n=list_nk_56n16n,list_NK_temp=list_NK_temp,list_nk_56hi=list_nk_56hi))
}


#----------------- Gating Baso---------------------------------------------------
gating_nk_56n16n_to_Baso<- function(){
  list_baso_poly<-lapply(1:length(gs),function(i){
    # adesso a partire dalla madre NK CD56-CD16- facciamo il gating ottenere le basophil cells.
    # prendiamo il flowFrame di nk5616n per i marker CD123 e HLA-DR
    
    temp <- fs.nk.56n16n[[i]]@exprs[, c(fluorochrome.chans['CD123'], fluorochrome.chans['HLA-DR'])]
    # originiariamente è stato usato nk.56n16n ovvero il cell population object  ma io uso il flowSet di nk.56n16n quindi cambio codice
    # se il numero di eventi di questo flowFrame è maggiore di 10
    if(nrow(fs.nk.56n16n[[i]])>10)
    {
      # flowPeaks = A fast and automatic clustering to classify the cells into subpopulations based 
      # on finding the peaks from the overall density function generated by K-means.
      
      # eseguiamo flowPeaks su temp per dividere la pop in clusters.
      # flowPeals.Res3 è un oggetto output di flowPeaks che contiene molti slot,vedi flow_cytometry_R_code per dettagli su flowPeaks.
      flowPeaks.Res3 <- flowPeaks(temp)
      # con print(flowPeaks.Res3) vediamo che flowPeaks ha individuato tre clusters(1,2,3),ogni cluster ha il suo id,ovvero 1,2,3.
      # visualizziamo temp in base a CD13: show(autoplot(temp,"PE-Cy7-A"))
      # picco grande intorno a 2,mini picchi a destra e sinistra.
      # calcoliamo il treshold per CD123 cercando di dividere la pop 50/50 (bimodal=T) e cerchiamolo nella coda della distribuzione.
      cd123.gate.nk<- deGate(fs.nk.56n16n[[i]],channel =  fluorochrome.chans['CD123'],bimodal=T,upper=T)
      # è pari a 2.964608,quindi dopo il picco grande
      # Ricorda:
      #  use.upper = T Logical.  If TRUE, forces to return the inflection point based on the first (last)
      # peak if upper=F (upper=T). default is F
      # upper = if TRUE, finds the change in the slope at the tail of the density curve, if FALSE,
      # finds it at the head.
      # Nota che la coda e la testa dipendono da dove si trova la media della distribuzione(valore più probabile della variabile),cioè se abbiamo un picco a destra ,la coda e' a sinistra e viceversa.
      # Quindi quando use.upper=T o =F estrare un preciso treshold nella coda(T) o nella testa(F)
      # Se use.upper=T e upper=F allora considera il primo picco della distribuzione
      # Se use.upper=T e upper=T considera l'ultimo picco della distribuzione 
      # Se use.upper=F torna il significato originale di upper
      
      #  prendiamo lo slot peaks dell'output di flowPeaks,di cui a sua volta prendiamo lo slot mu che contiene le media di tutte le cellule nei K clusters per ogni dimensione(in questa caso 2)
      # usiamo questa formula (flowPeaks.Res3$peaks$mu[, 1]  - 3)^2 + (flowPeaks.Res3$peaks$mu[, 2]  - 1)^2 e troviamo l'indice dell'elemento minore del risultato. In questo modo troviamo l'id del cluster delle basophil cells.
      # mu indica la media degli eventi per ogni k cluster(righe di mu) per ogni dimensione(colonne di mu)
      basophilcluster.id <- which.min((flowPeaks.Res3$peaks$mu[, 1]  - 3)^2 + (flowPeaks.Res3$peaks$mu[, 2]  - 1)^2)
      # prendiamo il flowFrame del CellPopulation object nk.56n16n
      basophils <- fs.nk.56n16n[[i]]
      # prendiamo la sua matrice di espressione e consideriamo lo slot cluster di peaks dove sono riportati i labels di ogni cell(a quale cluster appartiene ogni cellula)
      # cerchiamo quali cellule dello slot cluster appartengono all'id del cluster dei basofili ovvero le cellule che hanno il label del cluster dei basofili
      # con which otteniamo l'indice di queste cellule.
      basophils@exprs <-  basophils@exprs[which(flowPeaks.Res3$peaks.cluster %in% basophilcluster.id), ]
      # Calcoliamoi gates boundaries della basophils population come:
      # c(max(cd123.gate.nk,min(basophils@exprs[,c(fluorochrome.chans['CD123'])])), max(basophils@exprs[,c(fluorochrome.chans['HLA-DR'])]))
      # ovvero il treshld per cd123 è pari il massimo valore fra il threshold cd123.gate.nk(calcolato con flowDensity) e il minimo dei valori di espressione per CD123. In questo modo otteno il threshold più a destra fra i due(il mio scopo ricorda è di ottenere le cell CD123+)
      # il treshold per HLAD-DR corrisponde al massimo valore di espressione di HLA-DR.
      # questo perchè la matrice di espressione dell'oggetto basophil è già stata modificata in modo che contiene solo le cellule appartenenti al basophil cluster.
      # Quindi i confini del gate(quindi i threshold) sono gli estremi della  stessa popolazione sui sto applicando flowDensity(è solo per ottenere un valore più preciso di boundaries)
      # prendiamo dunque la popolazione a destra del treshold per CD123 e a sinistra di HLA_DR,quindi una pop CD123+Dr-
      basophils <- flowDensity(basophils, channels = c(fluorochrome.chans['CD123'], fluorochrome.chans['HLA-DR']), position = c(T, F), 
                               gates = c(max(cd123.gate.nk,min(basophils@exprs[,c(fluorochrome.chans['CD123'])])), max(basophils@exprs[,c(fluorochrome.chans['HLA-DR'])])), ellip.gate = T, scale = 0.975)
      # dunque riprendo la popolazione madre nk.56n16n,e applico flowDensity usando il filter rifinito che ho ottenut prima che si trova nel CellPopulation object basophils.
      basophils<- flowDensity(fs.nk.56n16n[[i]], channels = c(fluorochrome.chans['CD123'], fluorochrome.chans['HLA-DR']), position = c(T,T), filter = basophils@filter)
    }else{
      # se invece il numero di eventi è inferiore a 10,calcoliamo il threshold con deGate() prendendo nella coda  e pop 50/50.
      cd123.gate.nk<- deGate(fs.nk.56n16n[[i]], channel = fluorochrome.chans['CD123'],bimodal=T,upper=T)
      # prendiamo la pop a destra del treshold per CD123 e a sinistra del treshold di HLA-Dr calcolato come 0.5th percentile
      basophils <- flowDensity(fs.nk.56n16n[[i]], channels = c(fluorochrome.chans['CD123'], fluorochrome.chans['HLA-DR']), position = c(T,F),percentile=c(NA,.5)) 
    }
    # construct the basophils polygon object
    baso.poly <- polygonGate(filterId = "Baso",.gate  =basophils@filter)
    return(baso.poly=baso.poly)
  })
  # ------- add gates to the gating tree
  
  names(list_baso_poly)<-sampleNames(gs)
  nodeID0<-add(gs,list_baso_poly,parent="CD56-CD16-") 
  recompute(gs)
  return(list(gs=gs,list_baso_poly=list_baso_poly))
}


# cd123.gate<-averageGates(fsApply(fs.nk.neg,deGate,channels.ind[9],tinypeak.removal = 1/50,percentile=NA,upper=T),2)
# names(cd123.gate)<-sampleNames(fs.nk.neg)
# Baso<- fsApply(fs.nk.neg, function(x){
#   
#   baso <- flowDensity(x,channels.ind[c(9,1)],position = c(T,NA),gates=c(cd123.gate,NA))
#   plotDens(x,channels.ind[c(9,1)],cex=2)
#   lines(baso@filter,type="l",lwd=2)
#   return(baso)
# })

# ---------------------- Gating Granulocytes --------------------------------------
gating_gran_to_mgran_im_gran<-function(i){
  list_all_gran_poly <- lapply(1:length(gs),function(i){
    # adesso facciamo il gating a partire dai granulociti per ottenere le subpopulation dei granulociti:
    # mature granulocitye e immature granulocites.
    # Visualizziamo la pop dei granulociti in base a CD16: show(autoplot(fs.gran[[2]],"FITC-A")),show(autoplot(fs.gran[[1]],"BV786-A","FITC-A"))
    # grande picco a destra intorno a 3,quindi la coda è a sinistra dove c'è qualche mini-picco(molto molto piccoli)
    # calcoliamo il treshold di CD16 sulla pop generale dei granulociti(CD66+) che abbiamo trovato in gating precedenti.
    # lo posiziniamo come primo picco  della distribuzione(use.upper = T, upper = F)
    CD16granulo.gate <- deGate(fs.gran[[i]], channel = c(fluorochrome.chans['CD16']), use.upper = T, upper = F,alpha=0.9,tinypeak.removal = 0.8)
    # Nei PNG metto al posto di 0.90,metto 0.70
    CD16granulo.gate <- CD16granulo.gate * 0.80
    
    # è pari a 2.578116 quindi a sinistra del picco grande.
    #  prendiamo la pop a destra di questo treshold,quindi la pop CD16+
    temp.flowD <- flowDensity(fs.gran[[i]], channels = c(fluorochrome.chans['CD11b'], fluorochrome.chans['CD16']), position = c(NA, T), gates = c(NA, CD16granulo.gate))
    # visualizziamo la pop temp.flowD rispetto a Cd11b: show(autoplot(temp.flowD@flow.frame,"BV786-A"))
    # grande picco a destra intorno a 3,mini picchi nella coda di sinistra
    # calcoliamo il treshold su questa pop temp CD16+ come primo picco della distribuzione
    CD11bgranulo.gate <- deGate(temp.flowD, channel = c(fluorochrome.chans['CD11b']), use.upper = T, upper = F, alpha = 0.5,tinypeak.removal = 0.5)
    # Nei PNG files al posto di 2.4 metto 2.0
    if(CD11bgranulo.gate>=2.0){
      CD11bgranulo.gate <- deGate(temp.flowD, channel = c(fluorochrome.chans['CD11b']), use.upper = T, upper = F, alpha = 0.005,tinypeak.removal = 0.5)
      
    }
    # è pari a 1.914792 per il primo flowFrame quindi a sinistra del picco grande. 
    # Ripartiamo dalla popolazione generale dei granulociti prendiamo la popolazione a destra del treshol di CD11b e CD16,quindi otteniamo
    # la pop CD11b+Cd16+,ovvero i mature granulocytes o mature neutrophils
    maturegranulo <- flowDensity(fs.gran[[i]], channels = c(fluorochrome.chans['CD11b'], fluorochrome.chans['CD16']), position = c(T, T), gates = c(CD11bgranulo.gate, CD16granulo.gate))
    # dalla stessa madre prendiamo la pop CD11b-CD16- ,ovvero immature granulocites o immature neutrophils di tipo 1.
    immaturegranulo <- flowDensity(fs.gran[[i]], channels = c(fluorochrome.chans['CD11b'], fluorochrome.chans['CD16']), position = c(F, F), gates = c(CD11bgranulo.gate, CD16granulo.gate))
    # dalla stessa madre prendiamo anche la pop CD11b+CD16- ovvero i granulocites veri e propri senza i neutrophili
    CD11bposCD16neg.granulo<- flowDensity(fs.gran[[i]], channels = c(fluorochrome.chans['CD11b'], fluorochrome.chans['CD16']), position = c(T, F), gates = c(CD11bgranulo.gate, CD16granulo.gate))
    # dalla stessa madre prendiamo la pop CD11b-CD16+ il secondo tipo di immature neutrophils
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
                CD16granulo.gate=CD16granulo.gate))
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
  
  names(list_mature_gran_poly)<-sampleNames(gs)
  nodeID0<-add(gs,list_mature_gran_poly,parent="Granulocytes")
  names(list_immature_gran_poly)<-sampleNames(gs)
  nodeID0<-add(gs,list_immature_gran_poly,parent="Granulocytes")
  names(list_CD11bposCD16neg_gran_poly)<-sampleNames(gs)
  nodeID0<-add(gs,list_CD11bposCD16neg_gran_poly,parent="Granulocytes")
  names(list_CD11bnegCD16pos_gran_poly)<-sampleNames(gs)
  nodeID0<-add(gs,list_CD11bnegCD16pos_gran_poly,parent="Granulocytes")
  recompute(gs)
  fs.11p16p<-getData(gs,"CD11b+CD16+ Mature Neutrophils")
  return(list(fs.11p16p=fs.11p16p,gs=gs,list_CD16granulo_gate=list_CD16granulo_gate,list_CD11bgranulo_gate=list_CD11bgranulo_gate))
}


 


# ---------------------------- Gating CD64+ ----------------------------------------------------------------------
gating_mneutro_to_CD64pos<-function(i){
  list_all_CD64_poly <-lapply(1:length(gs),function(i){
    # infine dalla pop CD11b+CD16+,ovvero mature neutrophils  facciamo il gating per ottenere la pop CD64+
    # questo e' il gate che rym ha chiesto di aggiungere in piu' a mehrnoush,non e' presente nella gating strategy originale.
    
    # visualizziamo la pop in base a CD64: show(autoplot(fs.11p16p[[1]],"Alexa Fluor 700-A"))
    # grande picco vicino a 3,picchi piatti prima del grande picco,ovvero nella cosa
    # calcoliamo il treshold per CD64 di questa pop in modo da ottenere una pop 50/50 e posizionato nella coda
    cd64.gate<- deGate(fs.11p16p[[i]],fluorochrome.chans["CD64"],upper=F,bimodal = T,alpha = 0.8)
    # è pari a 1.682504 per il primo flowFrame,quindi prima del grande picco.
    # prendiamo la pop sinistra di questo threshold,quindi la pop CD64-,estriamo il filter dall CellPopulation object e facciamo con esso il PolygonGate.
    # --- construct polygon gate object
    CD64n.poly <- polygonGate(filterId ="CD64-",.gate = flowDensity(fs.11p16p[[1]],fluorochrome.chans[c("CD64","CD11b")],position = c(F,NA),gates=c(cd64.gate,NA))@filter)
    # stessa cosa per la pop a destra del treshold,quindi CD64+
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



#------ Data to return --------------------------------------------------------------------
generate_final_data<-function(){
  # scriviamo il codice per riportare creare l'output dei dati.
  # dalla librearia cytoUtils e flowWorkspace,usiamo la funzione getPopStats per ricavare info statistiche sul gatingHierarchy del gatingSet.
  # getPopStats returns a table population statistics for all populations in the gating hierarchy.
  # getProp returns the proportion of cells in the gate, relative to its parent
  # getTotal returns the total number of events included in this gate.
  
  counts<-getPopStats(gs, format="wide",statistic="freq")
  #print(counts)
  # usiamo la funzioneGatingSet2flowJo per creare un flowJo workspace 
  # GatingSet2flowJo = Convert a GatingSet to flowJo workspace 
  # in questo modo i dati possono essere visualizzati in flowJO.
  path.output_final_data <- paste0(path.output,"/","WSP","/")
  
  #GatingSet2flowJo(gs, outFile = paste0("/mnt/data/HIPC-EPIC/Myeloid/Gambia-0818/WSP/",sampleNames(fs)[1],"_AutomatedWSP.wsp"))
  GatingSet2flowJo(gs, outFile = paste0(path.output_final_data,sampleNames(fs)[1],"_AutomatedWSP.wsp"))
  # Ti ha prodotto due files di output nella cartella WSP uno di essi e' un file CSV che serve a flowJo quando aprira' il file .wsp per fare il gating dei clean,che indica i cluster di ogni cellula. 
  # se ti spunta l'errore "cannot find non boolean children under /Margin" o una cosa del genere,esso e' connesso al fatto quando hai aggiunto il nodo dei clean
  # tramite un vettore fattoriale hai dovuto indicare due nodi nell'argomento names (non ti trova i nodi in automatico il base al livello) perche' hai i pacchetti non aggiornati.
  # GatingSet2flowJo(gs, outFile = paste0("/mnt/data/HIPC-EPIC/Myeloid/Gambia-0818/WSP/",sampleNames(fs)[1],"_AutomatedWSP.wsp"))
  
}



#--------  Function to import and set up the dataset of manual counts
# pre_prop_manual_data <- function(){
# 
#   manual_data<- read.csv(file="/home/rstudio/data/Sample info/PNG_Pilot_wbc_counts_export_(all populations).csv",header=TRUE,sep=",")
# 
#   #------- riassegnamo nomi colonne del dataset manuale
#   # commentate sono le pop di cui ho il population object ma non le ho estratte dalla funzione
#   pop_names_manual <- colnames(manual_data)
#   pop_names_manual[which(pop_names_manual=="X")]<- "Sample"
#   pop_names_manual[which(pop_names_manual=="Size")]<- "All cells"
#   pop_names_manual[which(pop_names_manual=="Live")]<- "Live cells"
# 
#   pop_names_manual[which(pop_names_manual=="CD11b...CD16..Neut")]<- "CD11b-CD16+  Immature Neutrophils 2"
#   pop_names_manual[which(pop_names_manual=="CD11b...CD16..Neut.1")] <- "CD11b+CD16+ Mature Neutrophils"
#   pop_names_manual[which(pop_names_manual=="CD11b...CD16..Neut.2")] <- "CD11b+CD16- Granulocytes"
#   pop_names_manual[which(pop_names_manual=="CD11b...CD16..Neut.3")] <- "CD11b-CD16-  Immature Neutrophils 1"
#   pop_names_manual[which(pop_names_manual=="non.granulocytes")] <- "CD45+CD66- (Non Granulocytes)"
#   pop_names_manual[which(pop_names_manual=="DC...HLADR..CD14..")] <- "HLADR+ CD14-"
#   pop_names_manual[which(pop_names_manual=="B.cells")] <- "B cells"
#   pop_names_manual[which(pop_names_manual=="DC...CD14.HLADR..")] <- "HLADR- CD14-"
#   pop_names_manual[which(pop_names_manual=="CD3.")] <- "CD3-"
#   pop_names_manual[which(pop_names_manual=="CD3..CD56..CD16.")] <- "CD56-CD16-"
#   pop_names_manual[which(pop_names_manual=="Basophils")] <- "Baso"
#   pop_names_manual[which(pop_names_manual=="NK.CD56..CD16.")] <- "CD56-CD16+ NK"
#   pop_names_manual[which(pop_names_manual=="NK.CD56dim.CD16.")] <- "CD56 dim CD16+ NK-"
#   pop_names_manual[which(pop_names_manual=="NK.CD56dim.CD16..")] <- "CD56 dim CD16- NK"
#   pop_names_manual[which(pop_names_manual=="NK.CD56high.CD16.")] <- "CD56 Hi NK"
#   pop_names_manual[which(pop_names_manual=="T.cells")] <- "CD3+ T cells"
#   pop_names_manual[which(pop_names_manual=="gdT.cells")] <- "gd T cells"
#   pop_names_manual[which(pop_names_manual=="T.cells.CD16.CD56.")] <- "CD3+ T cells"
#   pop_names_manual[which(pop_names_manual=="T.cells.CD16.CD56..1")] <- "CD3+ T cells"
#   pop_names_manual[which(pop_names_manual=="T.cells..CD16..CD56.")] <- "CD3+ T cells"
#   pop_names_manual[which(pop_names_manual=="T.cells.CD16.CD56..2")] <- "CD3+ T cells"
#   pop_names_manual[which(pop_names_manual=="Monocytes")] <- "HLADR+ CD14+"
#   pop_names_manual[which(pop_names_manual=="Classical.Monocytes")] <- "Classical Monocytes"
#   #pop_names_manual[which(pop_names_manual=="Non.Classsical.monocytes")] <- ""
#   colnames(manual_data) <- pop_names_manual
#   #--- creiamo dataframe che associa a ogni nome di una pop la sua parente
#   list_child_vs_parent<-lapply(1:length(getNodes(gs)), function(i){
#     splitted_string <-str_split(getNodes(gs)[i],"/")
#     child_pop <- splitted_string[[1]][length(splitted_string[[1]])]
#     parent_pop <- splitted_string[[1]][length(splitted_string[[1]])-1]
#     return(list(child_pop=child_pop,parent_pop=parent_pop))
#   })
# 
#   list_all_pop_freq <- lapply(1:length(list_child_vs_parent), function(i){
#     child_pop_x <- list_child_vs_parent[[i]]$child_pop
#     if(child_pop_x %in% colnames(manual_data)){
#       parent_pop_x <- list_child_vs_parent[[i]]$parent_pop
#       if(parent_pop_x %in% colnames(manual_data)){
#         pop_x_freq <- manual_data[,child_pop_x]/manual_data[,parent_pop_x]
# 
#         return(list(pop_x_freq=pop_x_freq,name=child_pop_x))
#       }
#     }
#   })
# 
#   indx<-which(list_all_pop_freq=="NULL")
#   list_all_pop_freq<-list_all_pop_freq[-indx]
# 
#   names<-sapply(1:length(list_all_pop_freq),function(i){
#     return(list_all_pop_freq[[i]]$name)
#   })
#   list_all_pop_freq_only<-lapply(1:length(list_all_pop_freq),function(i){
#     return(list(pop=list_all_pop_freq[[i]]$pop_x_freq))
#   })
#   df <- as.data.frame(list_all_pop_freq_only)
#   colnames(df)<-names
#   manual_data_freq_on_parents <- df
#   samples<-manual_data$Sample
#   samples_split<-strsplit(as.character(samples),": ")
#   samples<-sapply(1:length(samples_split),function(i){
#     return(samples_split[[i]][2])
#   })
#   manual_data_freq_on_parents<-cbind(samples,manual_data_freq_on_parents)
#   return(manual_data_freq_on_parents)
# }














####################################### Execute functions ##################################
#Gating execution

main<-function(n_sample="All",path_comp_matrix,path_dir_files,path_fcs_files_fixing=NULL,plots=T,flowjo_wsp=F){
  #----- preprocessing
  # importiamo funzioni utili al pre-processing
  setwd("/home/rstudio/Code_Bcells_Myeloid_cells")
  path_utils <<- paste0("./Utils.R")
  source(path_utils)
  import_flowPrep()
  output_import_files<-import_files_v2(n_sample=n_sample,files_path=path_dir_files)

  fs<<-output_import_files$fs
  FCSfiles_path_list <<- output_import_files$FCSfiles_path_list
  path.output <<- create_path_output_dir("Myeloid_PNG")
  output_find_markers<<-find_markers_indices()
  scat.chans<<-output_find_markers$scat.chans
  fluorochrome.chans<<-output_find_markers$fluorochrome.chans
  Time.channel<<-output_find_markers$Time.channel

  gs<<- get_GatingSet()
  # --- gating dei marginal event
  output<<-margin_gating()
  fs.marg<<-output$fs_marg
  gs<<-output$gs
  #---- compensation e trasformation
  comp<<-import_comp_matrix_v2(path_comp_matrix=path_comp_matrix)
  out<<-comp_and_transform()
  gs<<-out$gs
  fs.marg<<-out$fs.marg
  # ----- compensation e trasformation
  output_prova<<-cleaning()
  fs.clean<<-output_prova$fs.clean
  gs<<-output_prova$gs
  clean.inds <<- output_prova$clean.inds
  #------ Singlets
  output_prova_2<<-gating_Time_to_singlets()
  fs.sngl<<-output_prova_2$fs.sngl
  gs<<-output_prova_2$gs
  singlets<<-output_prova_2$list_singlets
  #------- Size e beads
  output_prova_3<<-gating_singlets_to_Size_beads()
  fs.size<<-output_prova_3$fs.size
  gs<<-output_prova_3$gs
  beads_poly<<- output_prova_3$beads_poly_list
  size_poly <<- output_prova_3$size_poly_list
  beads <<-output_prova_3$list_beads

  # Rm("All cells",gs)
  # Rm("Beads",gs)

  #---- Live
  output_prova_4<<-gating_size_to_live()
  fs.live<<-output_prova_4$fs.live
  gs<<-output_prova_4$gs
  live_poly <<- output_prova_4$live_poly_list
  #plotGate(gs,"Live cells")
  #------ to Granulocytes e non gran
  output_prova_5<<-gating_live_to_granulocytes()
  fs.gran<<-output_prova_5$fs.gran
  fs.66n45p<<-output_prova_5$fs.66n45p
  gs<<-output_prova_5$gs
  gran_poly<<- output_prova_5$lista_poly_gran
  cd66_45_poly <<- output_prova_5$lista_poly_cd45
  # plotGate(gs,"Granulocytes")
  # plotGate(gs,"CD45+CD66- (Non Granulocytes)")
  # Rm("Granulocytes",gs)
  # Rm("CD45+CD66- (Non Granulocytes)",gs)
  #------ Gating to monocytes (14-) e 14+
  output_prova_6<<-gating_ngran_to_DRpos()
  fs.mono<<-output_prova_6$fs.mono
  fs.drn.14n<<-output_prova_6$fs.DRn.14n
  fs.drp.14n<<-output_prova_6$fs.DRn.14n
  gs<-output_prova_6$gs
  monocytes<<-output_prova_6$list_monocytes
  DRp_14n <<- output_prova_6$list_DRp_14n
  DRn_14n <<- output_prova_6$list_DRn_14n
  # Rm("HLADR+ CD14+",gs)
  # Rm("HLADR+ CD14-",gs)
  # Rm("HLADR- CD14-",gs)
  # plotGate(gs,"HLADR+ CD14+")
  # plotGate(gs,"HLADR+ CD14-")
  # plotGate(gs,"HLADR- CD14-")
  # ------ Gating to classical mono
  output_prova_7<<-gating_mono_to_classic_mono()
  gs<<-output_prova_7$gs
  fs.classmono <<- output_prova_7$fs.classmono
  CD16mon_gate <<- output_prova_7$list_CD16mon_gate

  #plotGate(gs,"Classical Monocytes")
  # Rm("Classical Monocytes",gs)
  #------- Gating to B_pdc_mdc
  output_prova_8<<-gating_drp_14n_to_B_pdc_mdc()
  gs<<-output_prova_8$gs
  # plotGate(gs,"B cells")
  # plotGate(gs,"mDC")
  # plotGate(gs,"pDC")
  # Rm("B cells",gs)
  # Rm("pDC",gs)
  # Rm("mDC",gs)

  mdc_poly <<- output_prova_8$lista_mdc_poly
  pdc_poly <<- output_prova_8$lista_pdc_poly
  bcell_poly <<- output_prova_8$lista_bcell_poly



  #------- Gating to t_gdt_cd3n
  output_prova_9<<-gating_drn_14n_to_T_gdt_cd3()
  fs.tcell<<-output_prova_9$fs.tcell
  fs.cd3n <<- output_prova_9$fs.cd3n
  gs<<-output_prova_9$gs
  # Rm("CD3+ T cells",gs)
  # Rm("CD3-",gs)
  # Rm("gd T cells",gs)
  # plotGate(gs,"CD3+ T cells")
  # plotGate(gs,"CD3-")
  # plotGate(gs,"gd T cells")

  CD3_gate <<- output_prova_9$list_CD3_gate
  gdt_poly <<- output_prova_9$list_gdt_poly
  #------- gating to NKT e NK
  output_prova_10<<-gating_nkt_nk()
  gs<<-output_prova_10$gs
  # plotGate(gs,"CD56 Hi NK")
  # plotGate(gs,"CD56-CD16-")
  # plotGate(gs,"CD56-CD16+ NK")
  # plotGate(gs,"CD56 dim CD16- NK")
  # plotGate(gs,"CD56 dim CD16+ NK")
  # Rm("CD56 Hi NK",gs)
  # Rm("CD56-CD16-",gs)
  # Rm("CD56-CD16+ NK",gs)
  # Rm("CD56 dim CD16- NK",gs)
  # Rm("CD56 dim CD16+ NK",gs)

  list_all_nkt_related_pops <<- output_prova_10$list_all_nkt_related_pops
  fs.nk.56n16n <<- output_prova_10$fs.nk.56n16n
  nk_56n16n <<- output_prova_10$list_nk_56n16n
  NK_temp <<- output_prova_10$list_NK_temp
  nk_56hi <<- output_prova_10$list_nk_56hi
  # ------ gating Baso
  output_prova_11 <<- gating_nk_56n16n_to_Baso()
  gs<<-output_prova_11$gs
  #plotGate(gs,"Baso")
  baso_poly <<- output_prova_11$list_baso_poly

  #---- gating  to mature e immature granulocytes
  output_prova_12 <<- gating_gran_to_mgran_im_gran()
  # Rm("CD11b+CD16+ Mature Neutrophils",gs)
  # Rm("CD11b-CD16-  Immature Neutrophils 1",gs)
  # Rm("CD11b+CD16- Granulocytes",gs)
  # Rm("CD11b-CD16+  Immature Neutrophils 2",gs)
  gs<<-output_prova_12$gs
  # plotGate(gs,"CD11b+CD16+ Mature Neutrophils")
  # plotGate(gs,"CD11b-CD16-  Immature Neutrophils 1")
  # plotGate(gs,"CD11b+CD16- Granulocytes")
  # plotGate(gs,"CD11b-CD16+  Immature Neutrophils 2")
  fs.11p16p<<-output_prova_12$fs.11p16p
  CD16granulo_gate <<- output_prova_12$list_CD16granulo_gate
  CD11bgranulo_gate <<- output_prova_12$list_CD11bgranulo_gate
  #------------ Gating to CD64+
  output_prova_13 <<- gating_mneutro_to_CD64pos()
  gs<<-output_prova_13$gs
  # plotGate(gs,"CD64+")
  # plotGate(gs,"CD64-")
  CD64_gate <<- output_prova_13$list_CD64_gate
  # Rm("CD64+",gs)
  # Rm("CD64-",gs)
  #---- Generate final Data
  if(flowjo_wsp==T){
    generate_final_data()
  }
  #---- Generates plots of each Gating step for each sample
  if(plots==T){
    generate_plots_myeloid()
  }
  return(gs)
}  

# #------------- section not usefull in final version: comparison with manual counts --------------------
# 
# #---- generate correlation plot
# manual_data_freq_on_parents<-pre_prop_manual_data()
# tab_stat_counts_freq<-pre_prop_automated_data()
# final_df<-make_final_dataframe()
# 
# corr_plot<-make_corr_plot(remove_p = "Baso|HLADR+ CD14+|HLADR- CD14-|CD56-CD16-|B cells|mDC|pDC")
# show(corr_plot)
# corr_plot<-make_corr_plot()
# 
# # fix: Baso,CD56-16- Nk,CD3-,CD3+ T cells,HLADR+,CD14+
# # ----- generate boxplots
# 
# boxplot <- generate_box_plots()

main("All",path_comp_matrix ="/home/rstudio/data/Comp Matrix_PNG_Feb2018-Updated.csv", 
     path_dir_files = "/home/rstudio/data/PNG pilot/")

