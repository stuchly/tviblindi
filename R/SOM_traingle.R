traingulation<-function(clusters,xdim,ydim) create_k_skeleton(TDA::alphaComplexFiltration(expand.grid(seq_len(xdim), seq_len(ydim))),clusters,2,type="rips")
