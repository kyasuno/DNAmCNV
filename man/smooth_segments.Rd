% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/smooth_segments.R
\name{smooth_segments}
\alias{smooth_segments}
\title{Smooth segments}
\usage{
smooth_segments(segs, loss.lrr = -0.1, gain.lrr = 0.1)
}
\arguments{
\item{segs}{list. Output of run_segmentation or annotate_segments.}

\item{loss.lrr}{numeric. Cutoff value of lrr for losses (default: -0.1). If NULL, segments are not classified.}

\item{gain.lrr}{numeric. Cutoff value of lrr for gains (default: 0.1). If NULL, segments are not classified.}
}
\value{
list. The same list with additional tbl named smoothed.
}
\description{
Merge consecutive segments in the same copy number state (loss, copy-neutral or gain).
}
