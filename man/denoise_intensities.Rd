% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/denoise_intensities.R
\name{denoise_intensities}
\alias{denoise_intensities}
\title{Denoise intensities using median filters recursively or non-recursively}
\usage{
denoise_intensities(x, chrArm, k, n.cores, recursive = FALSE, verbose = FALSE)
}
\arguments{
\item{x}{numeric vector.}

\item{chrArm}{a factor of chromosome arms of the same size as x.}

\item{k}{window size to be used for the sliding window (actually half-window size).}

\item{recursive}{logical. Perform median filters recursively or non-recursively (default: FALSE).}

\item{verbose}{logical.}
}
\value{
numeric vector.
}
\description{
Denoise intensities using median filters recursively or non-recursively
}
