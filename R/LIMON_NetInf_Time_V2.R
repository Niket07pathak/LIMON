

LIMON_NetInf_Time_V2 <- function(Obj, method = "glasso",
                              sel.criterion = "stars",
                              lambda.min.ratio,
                              nlambda,
                              pulsar.select = TRUE,
                              pulsar.params = list(),
                              icov.select = pulsar.select,
                              icov.select.params = pulsar.params,
                              lambda.log = TRUE) {
  
  library(SpiecEasi)
  
  CovMatrix_Time <- list()
  SpeicEasi_Time <- list()
  CovMatrix_Diff <- list()
  
  # Progress bar
  pb <- txtProgressBar(min = 0, max = length(Obj[["Corrected_Counts_Time"]]), style = 3)
  
  for (i in seq_along(Obj[["Corrected_Counts_Time"]])) {
    count_table <- Obj[["Corrected_Counts_Time"]][[i]]
    network_name <- paste0("Net_", i)
    cov_matrix_name <- paste0("Matrix_", i)
    
    if (method == "sparcc") {
      # Run SparCC directly
      net_data <- sparcc(count_table)
      matrix <- as.matrix(net_data$Cor)
      colnames(matrix) <- rownames(matrix) <- colnames(count_table)
    } else {
      # Run SPIEC-EASI
      net_data <- spiec.easi(data.matrix(count_table), method = method,
                             sel.criterion = sel.criterion, lambda.min.ratio = lambda.min.ratio,
                             nlambda = nlambda, pulsar.select = pulsar.select,
                             pulsar.params = pulsar.params, icov.select = icov.select,
                             icov.select.params = icov.select.params, lambda.log = lambda.log)
      
      # Extract optimal matrix
      if (method == "glasso") {
        matrix <- as.matrix(getOptCov(net_data))
      } else if (method == "mb") {
        matrix <- as.matrix(getOptBeta(net_data))
      } else {
        stop("Invalid method. Choose 'glasso', 'mb', or 'sparcc'.")
      }
      
      colnames(matrix) <- rownames(matrix) <- colnames(count_table)
    }
    
    # Store network and matrix
    SpeicEasi_Time[[network_name]] <- net_data
    CovMatrix_Time[[cov_matrix_name]] <- matrix
    
    # Update progress bar
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  
  # Store results
  Obj[["SpiecEasi_Time"]] <- SpeicEasi_Time
  Obj[["CovMatrix_Time"]] <- CovMatrix_Time
  
  # Pairwise difference computation
  N <- length(CovMatrix_Time)
  for (i in 1:(N - 1)) {
    for (j in (i + 1):N) {
      diff_table <- CovMatrix_Time[[paste0("Matrix_", j)]] - CovMatrix_Time[[paste0("Matrix_", i)]]
      diff_table_name <- paste0("Matrix_Diff_", j, i)
      CovMatrix_Diff[[diff_table_name]] <- diff_table
    }
  }
  Obj[["CovMatrix_Diff"]] <- CovMatrix_Diff
  Obj
}









# #testing
# L_obj3 <- LIMON_NetInf_Time_V2(Obj = L_obj2, method = "sparcc")

