.upsample.labels<-function(labels,N=10000,replace=TRUE,takeall="ungated"){
    out<-NULL
    out<-which(labels %in% takeall)
    for (i in unique(labels[!(labels %in% takeall)])) out<-c(out,sample(which(labels==i),size=ifelse(replace,N,min(N,length(which(labels==i)))),replace=replace))
    return(out[sample(length(out))])
}

##METHOD CHANGED - METHOD ADDED
merge_tviblindi<-function(x,fcsout="concatenated_fcs.fcs",normalize=NULL,selected_only=FALSE,shuffle=TRUE){
    stopifnot(is_list(x))
    labl<-unique(sapply(x,FUN=function(x) length(x$labels)))
    if (length(labl)>1) stop("Incompatible labeling")
    makefcs<-TRUE
    fcss<-NULL
    for (i in 1:length(x)){
        if (is.null(x[[i]]$fcs)) {
                               makefcs<-FALSE
                               break
        }
        fcss<-c(fcss,x[[i]]$fcs)
    }

    if(!makefcs){
        data<-NULL
        labels<-list()[1:(labl+1)]
        names(labels)<-c(names(x[[1]]$labels),"fileID")
        for (i in 1:length(x)){
            data<-rbind(data,x[[i]]$data)
            for (j in 1:labl) labels[[j]]<-c(labels[[j]],as.character(x[[i]]$labels[[j]]))
            labels[[labl+1]]<-c(labels[[labl+1]],rep(as.character(i),length(labels[[1]])))
        }
        x<-tviblindi(data=data,labels=labels)
        return(x)
    }

    data<-events_sel<-NULL
    labels<-list()[1:(labl+1)]
    Nf<-length(grep("fileID",names(x[[1]]$labels)))
    names(labels)<-c(names(x[[1]]$labels),paste("fileID",Nf,sep="_"))
    offset<-0
    fcs<-list()[1:length(x)]
    for (i in 1:length(x)){
        fcs[[i]]<-flowCore::read.FCS(fcss[i])
        data<-rbind(data,x[[i]]$data)
        for (j in 1:labl) labels[[j]]<-c(labels[[j]],as.character(x[[i]]$labels[[j]]))
        labels[[labl+1]]<-c(labels[[labl+1]],rep(as.character(i),length(x[[i]]$labels[[1]])))
        if (selected_only) ev_loc<-1:length(x[[i]]$events_sel) else ev_loc<-x[[i]]$events_sel
        events_sel<-c(events_sel,ev_loc+offset)
        offset<-offset+ifelse(selected_only,length(ev_loc),nrow(flowCore::exprs(fcs[[i]])))


    }

    if (!is.null(normalize)){
        if (normalize=="perc") data<-normalize.perc(data) else if (normalize=="scale") data<-scale(data)
    }
    if (shuffle) shuff<-sample(1:nrow(data)) else shuff<-1:nrow(data)
    if(!is.null(fcsout)){
        fcs<-.concat_fcs(fcs,x,params=paste("fileID",Nf,sep="_"),selected_only=selected_only)
        flowCore::write.FCS(fcs,filename=fcsout)
    }
    for (i in 1:(labl+1)) labels[[i]]<-labels[[i]][shuff]
    x<-tviblindi(data=data[shuff,],labels=labels,events_sel=events_sel[shuff],fcs=fcsout)
    return(x)

}

##adapted from package FlowCIPHE.
.concat_fcs<-function (flow.frames,x, params = "Flag",selected_only=FALSE)
{
    stopifnot(is_list(flow.frames))
    ff.concat <- NULL
    n <- length(flow.frames)
    for (i in 1:n) {
        if (selected_only) ff.raw <- flow.frames[[i]][x[[i]]$events_sel] else ff.raw <- flow.frames[[i]]
        p <- matrix(i, nrow = nrow(ff.raw), ncol = 1, dimnames = list(NULL,
            params))
        new.col <- as.vector(p)
        ff.raw <- .enrich_fcs(ff.raw, new.col, nw.names = params)
        if (is.null(ff.concat)) {
            ff.concat <- ff.raw
        }
        else {
            exprs(ff.concat) <- rbind(exprs(ff.concat), exprs(ff.raw))
        }
    }
    return(ff.concat)
}

.enrich_fcs<-function (original, new.column, nw.names = NULL)
{
    new_p <- flowCore::parameters(original)[1, ]
    new_p_number <- as.integer(dim(original)[2] + 1)
    rownames(new_p) <- c(paste0("$P", new_p_number))
    allPars <- BiocGenerics::combine(flowCore::parameters(original), new_p)
    if (is.null(nw.names)) {
        new_p_name <- "new_col"
    }
    else {
        new_p_name <- nw.names
    }
    allPars@data$name[new_p_number] <- new_p_name
    allPars@data$desc[new_p_number] <- new_p_name
    new_exprs <- cbind(original@exprs, new.column)
    colnames(new_exprs) <- c(colnames(original@exprs), new_p_name)
    new_kw <- original@description
    new_kw["$PAR"] <- as.character(new_p_number)
    new_kw[paste0("$P", as.character(new_p_number), "N")] <- new_p_name
    new_kw[paste0("$P", as.character(new_p_number), "S")] <- new_p_name
    new_kw[paste0("$P", as.character(new_p_number), "E")] <- "0,0"
    new_kw[paste0("$P", as.character(new_p_number), "G")] <- "1"
    new_kw[paste0("$P", as.character(new_p_number), "B")] <- new_kw["$P1B"]
    new_kw[paste0("$P", as.character(new_p_number), "R")] <- new_kw["$P1R"]
    new_kw[paste0("flowCore_$P", as.character(new_p_number),
        "Rmin")] <- new_kw["flowCore_$P1Rmin"]
    new_kw[paste0("flowCore_$P", as.character(new_p_number),
        "Rmax")] <- new_kw["flowCore_$P1Rmax"]
    new_fcs <- new("flowFrame", exprs = new_exprs, parameters = allPars,
        description = new_kw)
    return(new_fcs)
}
