#' Root phylogenetic tree of a phyloseq object
#'
#' @description This funtion defines the leaf with the longest path as the root of the phylogenetic tree.
#' This makes results reproducible by avoiding the behaviour of some functions that would otherwise pick a
#' random leaf as the root of an unrooted phylogenetic tree.
#' Based on answers in https://github.com/joey711/phyloseq/issues/597. The function requires the packages
#' 'ape' and 'data.table' to be installed.
#' @author Simeon Rossmann
#' @seealso Discussion and answers in [related GitHub thread](https://github.com/joey711/phyloseq/issues/597)
#'
#' @param physeq a phyloseq object containing a phylogenetic tree to be rooted in an outgroup.
#'
#' @return a rooted phylogenetic tree.
#'
#' @examples
#' phyloseq::phy_tree(ps) <- root_tree_in_outgroup(physeq = ps)
#'
#'@export
root_tree_in_outgroup <- function(physeq = ps){
  if(requireNamespace(c("ape", "data.table"), quietly = TRUE)){
    phylo_tree <- phyloseq::phy_tree(physeq)
    tips <- ape::Ntip(phylo_tree)
    tree_data <- base::cbind(
        data.table::data.table(phylo_tree$edge),
        data.table::data.table(length = phylo_tree$edge.length))[1:tips] %>%
      base::cbind(data.table::data.table(id = phylo_tree$tip.label))
    # longest terminal branch as outgroup
    out_group <- tree_data[which.max(length)]$id
    new_tree <- ape::root(phylo_tree, outgroup=out_group, resolve.root=TRUE)
    message("Tree successfully rooted.")
  }else{
    stop("The function 'root_tree_in_outgroup' requires the packages 'ape' and 'data.table' to be installed. Please make sure those packages can be loaded.")
  }
  return(new_tree)
}
