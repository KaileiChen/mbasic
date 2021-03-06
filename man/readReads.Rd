% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/InputDataProcessor.R
\name{readReads}
\alias{readReads}
\title{Read sequencing reads from either a BAM or a BED file.}
\usage{
readReads(reads, extended, fragLen = 200, pairedEnd = FALSE,
  use.names = FALSE, format)
}
\arguments{
\item{reads}{The sequencing file.}

\item{extended}{A boolean value for whether each read will be extended to the fragment length.}

\item{fragLen}{A numeric value for the fragment length.}

\item{pairedEnd}{A boolean value for whether the sequencing file is paired end.}

\item{use.names}{A boolean value to be passed to \link{GenomicRanges} functions.}

\item{format}{A string of file format. Must be either 'BAM' or 'BED'.}
}
\value{
A \linkS4class{RangedData} object.
}
\description{
Read sequencing reads from either a BAM or a BED file.
}
\author{
Samuel Younkin \email{syounkin@stat.wisc.edu}
}

