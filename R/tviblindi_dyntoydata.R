#' dyntoy datatset for tviblindi
#' 
#' dyntoy dataset is created as follows, principle component analisis is 
#' performed and first 25 PCs is used.
#' 
#' devtools::install_github('dynverse/dyntoy')
#' 
#' n_events <- 10000
#' n_features <- 3896
#' set.seed(12345)
#' d <- dyntoy::generate_dataset(
#' id           = 'tviblindi_dyntoy_test',
#'   model        = 'connected',
#'   num_features = n_features,
#'   num_cells    = n_events
#' )
#'
#'
#' @docType data
#'
#' @usage data(tviblindi_dyntoydata)
#'
#' @format An object of class \code{"data.frame"}; First column group_id 
#' character vector; columns 2:26 first 25 PCs
#'
#' @keywords datasets
#'
#'
#' @examples
#' data(tviblindi_dyntoydata)
#' group_is<-tviblindi_dyntoydata[,1]
#' datainput<-as.matrix(tviblindi_dyntoydata[,-1])
#' tv1<-tviblindi(data=data,labels=group_id)
#' DimRed(tv1)
#' DimRed(tv1,method="umap")
#'
#' Set_origin(tv1,label = "M4",origin_name = "M4_hitting_time")
#' Set_origin(tv1,label = "M4",origin_name = "M4_hitting_distance")

#' KNN(tv1)
#' Cluster(tv1) #kmeans clustering; K=625 clusters
#' Filtration(tv1) #default setting is too conservative, less simplices could be created with same resolution (e.g. Filtration(tv1,alpha2=1))
#' 
#' Pseudotime(tv1,weighted = FALSE,origin_name = "M4_hitting_time")
#' Walks(tv1,N=1000,origin_name = "M4_hitting_time")
#'
#'Pseudotime(tv1,weighted = TRUE,origin_name = "M4_hitting_distance")
#' Walks(tv1,N=1000,origin_name = "M4_hitting_distance")

#' launch_shiny(tv1)
#' times <- attr(grav, "time")
#' phe <- grav$pheno
"tviblindi_dyntoydata"