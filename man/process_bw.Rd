% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rate_estimation_helpers.R
\name{processBw}
\alias{processBw}
\title{Process bigwigs}
\usage{
processBw(bw, strand)
}
\arguments{
\item{bw}{\code{\link[GenomicRanges]{GRanges-class}} object representing
pro-seq counts}

\item{strand}{a string representing if counts are on the plus strand with
+' or the minus
strand with '-'}
}
\value{
A \code{\link[GenomicRanges]{GenomicRanges-class}} object with
single basepair
resolution
}
\description{
Creates GRanges object with strand information and counts set to single
basepair resolution and
absolute value. Keeps only standard chromosome information and excludes
count regions set to 0.
}
