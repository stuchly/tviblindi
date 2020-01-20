sample_points<-function(ws,  nbL = 100L) {
    out<-.Call('_tviblindi_internal_sample_points', PACKAGE = 'tviblindi', ws, nbL)
    out<-do.call("rbind",out$landmarks)
    out
}
