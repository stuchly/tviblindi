tviblindi_from_flowjo<-function(fcsn,wsp,fcstable,eventsufix="_G[0-9]+$",eventprefix="/[A-Z]_",origin="^A_",normalize=FALSE,sep="\t"){
  require(CytoML)
  require(flowWorkspace)

  ws<-open_flowjo_xml(wsp)
  protilatky<-read.table(fcstable,header=TRUE,sep=sep)
  fcs<-flowCore::read.FCS(fcsn)
  chans<-protilatky$name[protilatky$use==1]
  chans<-which(fcs@parameters@data$name %in% chans)
  gs<-flowjo_to_gatingset(ws,name=1,subset=1)
  PP<-gs_get_pop_paths(gs[[1]])
  ppi<-grep(eventsufix,PP)
  cat("events used: ", PP[ppi],"\n")
  events_sel<-NULL;for (i in ppi) events_sel<-c(events_sel,which(gh_pop_get_indices(gs[[1]],PP[i[1]])))
  COO<-exprs(fcs)[events_sel,chans]
  colnames(COO)<-as.character(fcs@parameters@data$desc[chans])
  COO<-asinh(COO/5.0)
  if (normalize=="perc") COO<-normalize.perc(COO) else if (normalize=="scale") COO<-scale(COO)
  ##label events
  leafs<-NULL
  for (ii in 1:length(PP)) if (length(grep(PP[ii],PP,fixed=TRUE))==1) leafs<-c(leafs,ii)
  leafs<-leafs[-1]
  leafs<-leafs[grep(eventprefix,PP[leafs])]
  indM<-NULL

  for (ii in leafs) indM<-cbind(indM,as.numeric(gh_pop_get_indices_mat(gs[[1]],PP[ii])))

  colnames(indM)<-gsub(".*/","",PP[leafs])
  gate_ev<-colnames(indM)[apply(indM,MARGIN=1,which.max)]
  ss<-which(rowSums(indM)==0)
  gate_ev[ss]<-"ungated" ##negative gate

  gate_ev_non<-gate_ev[events_sel]
  gate_ev_non<-as.factor(gate_ev_non)
  stems           <- levels(gate_ev_non)[grep(origin,levels(gate_ev_non))]

  tv1<-tviblindi(COO,labels = gate_ev_non,fcs_path = fcsn,events_sel = events_sel)
  Set_origin(tv1,stems)
  if(is.null(tv1$origin)) stop("Wrong population of origin")
  print(tv1)
  return(tv1)
}
