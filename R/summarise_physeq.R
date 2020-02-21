#' Summarize a phyloseq object
#'
#' @author Simeon Rossmann
#' @param physeq a phyloseq object to be summarized
#' @param ASV_sublist a character vector of ASVs/OTUs (phyloseq taxa). Empty by default.
#' If provided the summary will contain the percentage of abundance the sublist comprises.
#' Compatible with cuphyr::list_subset_ASVs().
#' @param sublist_id a character string describing the subset of taxa provided in ASV_sublist
#' @param samp_names TRUE/FALSE. If TRUE (default), sample names in the subset are included in the summary.
#'
#' @description Summarize the total number of ASVs, number of ASVs in a subset and percentage
#' of reads this ASV subset contains (of total) of a given phyloseq object.
#' @return a character string summarising the provided phyloseq object
#'
#' @examples summarise_physeq(ps)
#' summarise_physeq(ps, ASV_sublist= c("ASV1"), sublist_id ="ASV1", samp_names=FALSE)
#'
#' phyt_ASVs <- cuphyr::list_subset_ASVs(physeq = ps, subv = c("Phytophthora","Phytopythium"), taxlvlsub="Genus")
#' summarise_physeq(ps, ASV_sublist= phyt_ASVs, sublist_id ="Phytophthora and Phytopythium", samp_names=FALSE)
#'
#' @export
summarise_physeq <- function(physeq = ps, ASV_sublist= c(), sublist_id = "sublist_id_not_defined", samp_names=TRUE){
  samps <- paste(phyloseq::sample_names(physeq), collapse = "\n")
  numTotalASVs <- length(phyloseq::taxa_names(physeq))
  obj <- deparse(substitute(physeq))

  # If a sublist is given, get the percentage of
  if(!is.null(ASV_sublist)){
    numberSubASVs<- length(ASV_sublist)
    ps_prcnt <- phyloseq::transform_sample_counts(physeq, function(ASV) (ASV/sum(ASV))*100)
    ps_prcnt_sub <- suppressWarnings(phyloseq::prune_taxa(ASV_sublist, ps_prcnt))
    samp_subPercent <- mean(phyloseq::sample_sums(ps_prcnt_sub))
    samp_subPercent <- format(round(samp_subPercent, 2), nsmall = 2)
  }

  #Compose substrings depending on arguments provided to the function
  strone <- paste("There are", numTotalASVs,"ASVs in the phyloseq object", paste0("'", obj,"'."))
  if(!is.null(ASV_sublist)){
    strtwo <- paste("Of this,", numberSubASVs, ifelse(numberSubASVs==1, "belongs", "belong"),
                    "to the provided subset", paste0("(",sublist_id,"),"), "representing",
                    samp_subPercent, "percent of abundance per sample on average.")
  }else{
    strtwo <- c()
  }
  if(samp_names){
    strthree <- paste("The samples in", paste0("'", obj,"'"), "are:\n", samps)
  }else{
    strthree <- c()
  }

  #Compose final output
  summaryPhyseq <- paste(strone, strtwo, strthree, sep = "\n")
  return(cat(summaryPhyseq))
}
