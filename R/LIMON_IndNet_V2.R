


LIMON_IndNet_V2 <- function(Obj, method = "glasso",
                            sel.criterion = "stars",
                            lambda.min.ratio,
                            nlambda,
                            pulsar.select = TRUE,
                            pulsar.params = list(),
                            icov.select = pulsar.select,
                            icov.select.params = pulsar.params,
                            lambda.log = TRUE) {
  
  library(SpiecEasi)
  
  Individual_Networks <- list()
  
  pb <- txtProgressBar(min = 0, max = length(Obj[["Corrected_Counts_Time"]]), style = 3)
  
  for (i in seq_along(Obj[["Corrected_Counts_Time"]])) {
    count_table <- Obj[["Corrected_Counts_Time"]][[i]]
    nsamples <- nrow(count_table)
    
    if (method == "sparcc") {
      net_all <- sparcc(count_table)$Cor
    } else {
      net_all <- as.matrix(getOptCov(Obj[["SpiecEasi_Time"]][[i]]))
    }
    
    for (j in 1:nsamples) {
      samplename <- rownames(count_table)[j]
      count_minus_j <- count_table[-j, , drop = FALSE]
      
      if (method == "sparcc") {
        net_minus_j <- sparcc(count_minus_j)$Cor
        net_individual <- nsamples * (net_all - net_minus_j) + net_minus_j
      } else {
        net_single <- spiec.easi(data = as.matrix(count_minus_j), method = method,
                                 sel.criterion = sel.criterion, lambda.min.ratio = lambda.min.ratio,
                                 nlambda = nlambda, pulsar.select = pulsar.select,
                                 pulsar.params = pulsar.params, icov.select = icov.select,
                                 icov.select.params = icov.select.params, lambda.log = lambda.log)
        
        net_minus_j <- as.matrix(getOptCov(net_single))
        net_individual <- nsamples * (net_all - net_minus_j) + net_minus_j
      }
      
      colnames(net_individual) <- rownames(net_individual) <- colnames(count_table)
      Individual_Networks[[samplename]] <- net_individual
    }
    
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  Obj[["Individual_Networks"]] <- Individual_Networks
  Obj
}

# Obj <- LIMON_IndNet_V2(Obj = L_obj3, method = "sparcc")


# # Extract edges and centralities
# L_obj7<- LIMON_IndEdges(Obj, threshold = 0.2)
# L_obj8 <- LIMON_Centralities(L_obj7, threshold = 0.2)
