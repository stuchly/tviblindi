suppressPackageStartupMessages(library(tviblindi))

data(tviblindi_dyntoydata)
group_id <- tviblindi_dyntoydata[, 1]
datainput <- as.matrix(tviblindi_dyntoydata[, -1])

tv1 <- tviblindi(data = datainput, labels = group_id)
Set_origin(tv1, label = "M4", origin_name = "M4_hitting_time")
KNN(tv1)
Cluster(tv1, K = 225)
Filtration(tv1)
Pseudotime(tv1, weighted = FALSE, origin_name = "M4_hitting_time")
Walks(tv1, N = 1000, origin_name = "M4_hitting_time")
DimRed(tv1, method = "umap")

cat("=== tv1 built, launching shiny on 0.0.0.0:3838 ===\n")

options(shiny.host = '0.0.0.0', shiny.port = 3838)
launch_shiny(tv1)
