#' Phyloseq to LIMON
#'
#' Create a formatted objects to run LIMON from a Phyloseq object
#'
#' @param ps_object a phyloseq object
#'
#' @returns Returns a LIMON Object with a Count and SampleData table in the correct format
#'
#' @export
#'
phyloseq_to_LIMON_Obj <- function(ps_object){

  # Exctract Counts
  Counts_ps <- as.data.frame(ps_object@otu_table)
  SampleData_ps <- as.data.frame(ps_object@sam_data)
  # Check if row names are identical and in the same order in Counts and SampleData
  if (!identical(rownames(Counts_ps), rownames(SampleData_ps))) {
    Counts_ps <- as.data.frame(t(Counts_ps))
  }

  list(Counts = Counts_ps,
       SampleData =SampleData_ps)
}
