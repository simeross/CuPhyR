#' Returns a vector of ASVs belonging to groups at a taxonomic level
#' @description A function to list ASVs belonging to a set of taxonomic terms
#' @author Simeon Rossmann
#' @param physeq a phyloseq object to extract ASVs/OTUs ('tax_names') from
#' @param subv a character vector of taxonomic terms at the level of interest
#' @param taxlvlsub a character string containing a taxonomic level present in the database
#'
#' @examples phyt_ASVs <- list_subset_ASVs(physeq = phyloseq_oomycetes, subv = c("Phytophthora","Phytopythium"), taxlvlsub="Genus")
#' @return a character vector of ASVs
#' @export
list_subset_ASVs <- function(physeq=ps, subv=c("e"), taxlvlsub="Kingdom"){
  #Preparing lists of the chosen taxonomic level and searching for overlap with values in vector
  taxlist <- phyloseq::tax_table(physeq)[,taxlvlsub]
  taxcounts <- stats::aggregate(data.frame(count=taxlist), list(value=taxlist), length)
  unilist <- taxcounts$value
  hitlist <- unilist[unilist %in% grep(paste0(subv, collapse = "|"), unilist, value = T)]
  #Searching for ASVs that contain 'subv' values at 'taxlvlsub' and building subset
  ASVlist <- character(0)
  for (taxgroup in hitlist){
    ASVlist <- c(ASVlist, rownames(subset(phyloseq::tax_table(physeq), phyloseq::tax_table(physeq)[, taxlvlsub] == taxgroup)))
  }
  return(ASVlist)
}
