#' Individual Network Inferences with Sparcc and Precomputed Lambda Support
#'
#' Function to perform network inferences following the method in lionessR. 
#' Supports SpiecEasi ("glasso", "mb") and SparCC. Also supports reusing 
#' precomputed lambda values per time point to speed up individual inference.
#'
#' @param Obj LIMON Object from LIMON_NetInf_Time()
#' @param method Network inference method: "glasso" [default], "mb", or "sparcc"
#' @param sel.criterion Method for model selection (SpiecEasi only). Either 'stars' [default] or 'bstars'
#' @param lambda.min.ratio Scaling factor for minimum lambda (SpiecEasi only)
#' @param nlambda Number of lambda values (SpiecEasi only)
#' @param pulsar.select Default TRUE (SpiecEasi only)
#' @param pulsar.params Default empty list()
#' @param icov.select Default "pulsar.select"
#' @param icov.select.params Default "pulsar.params"
#' @param lambda.log Default TRUE
#' @param lambda.per.time Optional list of precomputed optimal lambdas (one per timepoint)
#'
#' @return LIMON object with individual networks per sample stored in Individual_Networks by sample name
#' @export
#'
LIMON_IndNet_V2 <- function(Obj, method = "glasso",
                            sel.criterion = "stars",
                            lambda.min.ratio,
                            nlambda,
                            pulsar.select = TRUE,
                            pulsar.params = list(),
                            icov.select = pulsar.select,
                            icov.select.params = pulsar.params,
                            lambda.log = TRUE,
                            lambda.per.time = NULL) {

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
        # Check for precomputed lambda
        if (!is.null(lambda.per.time)) {
          net_single <- spiec.easi(data = as.matrix(count_minus_j),
                                   method = method,
                                   lambda = lambda.per.time[[i]],
                                   pulsar.select = FALSE)
        } else {
          net_single <- spiec.easi(data = as.matrix(count_minus_j),
                                   method = method,
                                   sel.criterion = sel.criterion,
                                   lambda.min.ratio = lambda.min.ratio,
                                   nlambda = nlambda,
                                   pulsar.select = pulsar.select,
                                   pulsar.params = pulsar.params,
                                   icov.select = icov.select,
                                   icov.select.params = icov.select.params,
                                   lambda.log = lambda.log)
        }

        if (method == "glasso") {
          net_minus_j <- as.matrix(getOptCov(net_single))
        } else if (method == "mb") {
          net_minus_j <- as.matrix(getOptBeta(net_single))
        }

        net_individual <- nsamples * (net_all - net_minus_j) + net_minus_j
      }

      colnames(net_individual) <- rownames(net_individual) <- colnames(count_table)
      net_name <- paste0(samplename, "_Time", i)
      Individual_Networks[[net_name]] <- net_individual
    }

    setTxtProgressBar(pb, i)
  }

  close(pb)
  Obj[["Individual_Networks"]] <- Individual_Networks
  Obj
}



# # Extract lambdas from LIMON_NetInf_Time object
# lambda.per.time <- lapply(L_obj3[["SpiecEasi_Time"]], function(x) {
#   x$lambda[x$select$stars$opt.index]
# })

# # Call updated function
# L_obj3 <- LIMON_IndNet_V2(Obj = L_obj3,
#                           method = "glasso",
#                           lambda.per.time = lambda.per.time,
#                           pulsar.select = FALSE)  # skip tuning
