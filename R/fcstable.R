fcstable<-function(fcs,table_file="channels.txt"){
    fcs<-read.FCS(fcs)
    res<-fcs@parameters@data
    res<-data.frame(res,use=0)
    write.table(res,file=table_file,col.names=TRUE,row.names=FALSE,sep="\t")
    return(invisible(res))
}
