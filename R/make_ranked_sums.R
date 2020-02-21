#' Rank ASV abundance sums for a given phyloseq object
#'
#' @description A function that ranks all samples in a phyloseq object by the sum of reads assigned to all ASVs in that sample (Abundance sum/Total ASV counts).
#' @author Simeon Rossmann
#' @param physeq a phyloseq object to rank by abundance
#' @param myset optional: a character string describing the ranked physeq object.
#'
#' @return a tibble that contains all sample metadata (fcts), their abundance (dbl),
#' their ranks (int) by sum of abundance, their original order (int) and a custom descriptor (chr).
#'
#' @examples ranked_controls <- make_ranked_sums(physeq = ps_contrls)
#'
#'@export
make_ranked_sums <- function(physeq = ps, myset = "Not specified") {
  #Parse metadata
  metadat <- phyloseq::sample_data(physeq)
  metadat_tbl <- suppressWarnings(metadat %>%
                                    dplyr::as_tibble()) %>%
    dplyr::mutate(Sample = rownames(metadat)) %>%
    dplyr::select(Sample, dplyr::everything())
  #Parse Sample counts
  totals_samps <- phyloseq::sample_sums(physeq)
  names_samps <- phyloseq::sample_names(physeq)
  sample_tbl <- dplyr::tibble(Sample = names_samps, Abundance = totals_samps[Sample])
  #Combine and rank
  ranked <- dplyr::left_join(metadat_tbl, sample_tbl, by="Sample") %>%
    dplyr::mutate(Order = 1:dplyr::n()) %>%
    dplyr::arrange(dplyr::desc(Abundance)) %>%
    dplyr::mutate(Rank = 1:dplyr::n()) %>%
    dplyr::arrange(Order) %>%
    dplyr::mutate(Set = myset)
  return(ranked)
}
