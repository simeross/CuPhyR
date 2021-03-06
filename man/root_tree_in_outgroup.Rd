% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/root_tree_in_outgroup.R
\name{root_tree_in_outgroup}
\alias{root_tree_in_outgroup}
\title{Root phylogenetic tree of a phyloseq object}
\usage{
root_tree_in_outgroup(physeq = ps)
}
\arguments{
\item{physeq}{a phyloseq object containing a phylogenetic tree to be rooted in an outgroup.}
}
\value{
a rooted phylogenetic tree.
}
\description{
This funtion defines the leaf with the longest path as the root of the phylogenetic tree.
This makes results reproducible by avoiding the behaviour of some functions that would otherwise pick a
random leaf as the root of an unrooted phylogenetic tree.
Based on answers in https://github.com/joey711/phyloseq/issues/597. The function requires the packages
'ape' and 'data.table' to be installed.
}
\examples{
phyloseq::phy_tree(ps) <- root_tree_in_outgroup(physeq = ps)

}
\seealso{
Discussion and answers in [related GitHub thread](https://github.com/joey711/phyloseq/issues/597)
}
\author{
Simeon Rossmann
}
