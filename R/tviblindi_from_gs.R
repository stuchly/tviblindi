require(CytoML)
require(flowWorkspace)

# wsp="20221125 BM MNCs DRR day 0.wsp"
# ws <- open_flowjo_xml(wsp)
# gs <- flowjo_to_gatingset(ws,name=1)
# for (i in 1:length(gs)) print(keyword(gs[[i]], "Category"))

tviblindi_from_gatingset<-function(gs,eventsufix = "_G$", eventprefix = "/[A-Z][0-9]*_", 
                                   origin = "^A_", normalize = FALSE,fcs_paths=NULL,merge_fcs=NULL){
    
    data<-NULL
    events_sel_all<-NULL
    gate_all<-list()[1:(length(eventprefix)+1)]
    Nevents<-0
    if (is.null(fcs_paths)) fcss<-gsub("_[0-9]+$","",sampleNames(gs)) else fcss<-fcs_paths
    for (iframe in 1:length(gs)){
        #print(fcss[iframe])
        PP <- gs_get_pop_paths(gs[[iframe]])
       
        ppi <- grep(eventsufix, PP)
        if (length(ppi) == 0) {
            warning(paste("Suffix not found in",fcss[iframe],"- returning NULL", "\n"))
            next
        }
        cat("events used: ", PP[ppi], "\n")
        events_sel <- NULL
        for (i in ppi) events_sel <- c(events_sel, which(gh_pop_get_indices(gs[[iframe]], 
                                                                            PP[i[1]])))
        
        events_sel <- unique(events_sel)
       
        efcs<-flowCore::exprs(gh_pop_get_data(gs[[iframe]]))
        
        if (iframe>1) events_sel_all<-c(events_sel_all,events_sel+Nevents) else events_sel_all<-events_sel
        Nevents<-nrow(efcs)+Nevents
        print(length(events_sel))
        efcs<-efcs[events_sel,]
        data<-rbind(data,efcs)
        gate_ev_non <- list()[1:length(eventprefix)]
        j <- 0
        for (pref in eventprefix) {
            j <- j + 1
            if (pref=="*"){
                leafs<-NULL
                for (i in 1:length(PP)) {
                    if (length(grep(PP[i],PP,fixed=TRUE))==1) leafs<-c(leafs,i)
                    
                }
                leafs<-leafs[-1]
                
            } else {
                leafs <- grep(pref, PP) 
            }
            indM <- NULL
            for (ii in leafs) indM <- cbind(indM, as.numeric(gh_pop_get_indices_mat(gs[[iframe]], 
                                                                                    PP[ii])))
            colnames(indM) <- gsub(".*/", "", PP[leafs])
            
            gate_ev <- colnames(indM)[apply(indM, MARGIN = 1, which.max)]
            ss <- which(rowSums(indM) == 0)
            gate_ev[ss] <- "ungated"
            gate_ev_non[[j]] <- gate_ev[events_sel]
        }
        gate_ev_non[[j+1]]<-rep(iframe,length(gate_ev_non[[j]]))
        for (j in 1:length(gate_ev_non)) gate_all[[j]]<-as.character(c(gate_all[[j]],gate_ev_non[[j]]))
    }
    
    #fcss<-paste(dir()[c(3,4,5,10)],gsub("_[0-9]+$","",sampleNames(gs)),sep="/")
    if (!is.null(merge_fcs) & length(fcss>1)){
      fcs<-list()[1:length(fcss)]
      for (i in 1:length(fcss)) fcs[[i]]<-flowCore::read.FCS(fcss[i])
      
      fcs <- tviblindi:::.concat_fcs(fcs, x, params = paste("fileID", 0, 
                                                sep = "_"), selected_only = FALSE)
      flowCore::write.FCS(fcs, filename = merge_fcs)
      
    }
    if (length(fcss)==1) fcsout=fcss
    
    names(gate_all)<-paste("labels",1:length(gate_all))
    #return(list(data,gate_all,events_sel_all))
    tv1<-tviblindi(data=data,labels=gate_all,events_sel = events_sel_all,fcs_path = merge_fcs)
    return(tv1)
    
    
}

