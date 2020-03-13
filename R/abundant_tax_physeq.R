#' Return top n most abundand taxonomic terms at a given level
#'
#' @author Simeon Rossmann
#' @description This function extracts the top n most abundant taxonomic terms at a given level.
#' It uses dplyr to calculate the most abundant elements.
#' Output_format will take "ps" or "tops" and either return the phyloseq object or a character vector of the n most common terms.
#' Ignores NA's by default but can show them (see parameters).
#' @param physeq a phyloseq object to rank
#' @param lvl a character vector defining the taxonomic level to rank
#' @param top an integer n defining the number of returned elements (e.g. 'top=5')
#' @param output_format "tops" or "ps" (default), depending on this value, function returns either a character vector or a pruned phyloseq object.
#' @param ignore_na TRUE/FALSE, default FALSE. If TRUE, excludes ASVs without assigned taxonomic term in the output.
#' @param silent TRUE/FALSE, default TRUE. If FALSE, prints the top n taxonomix terms as a message.
#'
#'
#' @examples ps_top10genera <- abundant_tax_physeq(physeq = ps, lvl="Genus", top=10)
#' top5species <- abundant_tax_physeq(physeq = ps, lvl="Species", top=5, output_format="tops", ignore_na = FALSE)
#'
#' @return a phyloseq object or character vector, depending on value of "output_format"
#' @export
abundant_tax_physeq <- function(physeq = ps, lvl="Class", top=10, output_format="ps", ignore_na = FALSE, silent=TRUE){
  #This function
  if(!exists("ps_tbl")){
    ps_tbl <- phyloseq::psmelt(physeq)
    warning("The tbl version of the phyloseq object was not found. It will be generated, which may take a while.
            To significantly increase the speed of this function, make a tbl of your physeq object
            and call it 'ps_tbl' using \n \n
            ps_tbl <- dplyr::as_tibble(phyloseq::psmelt(ps))")
  }
  #Some quotation magic to get around the evaluation rules in the tidyverse (first unquote 'lvl', then enquo() for later evaluation with !!)
  lvlx <- dplyr::sym(lvl)
  lvlx <- dplyr::enquo(lvlx)
  #Group at lvl, sort abundance sums per group descending, then extract the top n and pull the taxonomic terms to vector.
  if(ignore_na){
    top_x <-  suppressWarnings(ps_tbl %>%
                                 dplyr::group_by(Tax = !!lvlx) %>%
                                 dplyr::summarise(Abundance=sum(Abundance))  %>%
                                 dplyr::arrange(dplyr::desc(Abundance)) %>%
                                 dplyr::filter(!is.na(Tax)) %>%
                                 dplyr::top_n(top) %>%
                                 dplyr::pull(Tax) %>%
                                 as.character())
  }else{
    top_x <-  suppressWarnings(ps_tbl %>%
                                 dplyr::group_by(Tax = !!lvlx) %>%
                                 dplyr::summarise(Abundance=sum(Abundance))  %>%
                                 dplyr::arrange(desc(Abundance)) %>%
                                 dplyr::top_n(top) %>%
                                 dplyr::pull(Tax) %>%
                                 as.character())
  }
  if(!silent){
    message(cat("\nThe top ", top, " most abundant annotated groups at the taxonomic level '", lvl,
                "' are:\n", sep =""), cat(top_x, sep="\n"))}
  #Return based on output_format
  if(output_format=="ps"){
    if(nrow(phyloseq::tax_table(physeq)[,lvl]) <= top){
      ps_topnt <- physeq
    }else{
      #Loop through all elements in top_x and extract the ASVs corresponding to the top n taxgroups at lvl from the specified physeq object. Then prune physeq.
      toptaxASVs <- suppressWarnings(ps_tbl %>%
                    dplyr::group_by(OTU, !!lvlx) %>%
                    dplyr::summarise() %>%
                    dplyr::filter(!!lvlx %in% top_x) %>%
                    dplyr::select(OTU) %>%
                    unlist() %>%
                    as.character())
      ps_topnt <- phyloseq::prune_taxa(toptaxASVs, physeq)}
    return(ps_topnt)}
  if(output_format=="tops"){
    return(top_x)

  }else{
    warning("Choose \"ps\" or \"tops\" for output_format! ", output_format, " is not a legitimate value.
            Phyloseq object returned as default")
    if(nrow(phyloseq::tax_table(physeq)[,lvl]) <= top){
      ps_topnt <- physeq
    }else{
      #Loop through all elements in top_x and extract the ASVs corresponding to the top n taxgroups at lvl from the specified physeq object. Then prune physeq.
      toptaxASVs <- suppressWarnings(ps_tbl %>%
                                       dplyr::group_by(OTU, !!lvlx) %>%
                                       dplyr::summarise() %>%
                                       dplyr::filter(!!lvlx %in% top_x) %>%
                                       dplyr::select(OTU) %>%
                                       unlist() %>%
                                       as.character())
      ps_topnt <- phyloseq::prune_taxa(toptaxASVs, physeq)}
    return(ps_topnt)}
}
