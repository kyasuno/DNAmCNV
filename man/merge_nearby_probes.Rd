% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/merge_nearby_probes.R
\name{merge_nearby_probes}
\alias{merge_nearby_probes}
\title{Merge nearby probes}
\usage{
merge_nearby_probes(crs, min.gapwidth, max.probes, n.cores = 1L)
}
\arguments{
\item{crs}{list. Output of \code{tangent_normalization}.}

\item{min.gapwidth}{integer. For \code{GenomicRanges::reduce} function (default: 500L).}

\item{max.probes}{integer. Maximum number of probes to be merged within a set of probes. Probes are split into nearly equally
\code{ceiling(k/max.probes)} groups, where k is the number of probes in the group.}
}
\value{
tbl. seqnames, start, end, width, nProbes, chrArm and median log R ratio (lrr).
}
\description{
Merge nearby probes
}
