% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/proteinSetReport.R
\name{runPianoGSEAnalysis}
\alias{runPianoGSEAnalysis}
\title{Run piano GSEA analysis}
\usage{
runPianoGSEAnalysis(geneProtIdent, pvalues, foldChanges = NULL, pianoGSC,
  minGSSize = 5, maxGSSIze = 300)
}
\arguments{
\item{geneProtIdent}{A list of proteins or gene identifiers, which matches the given gene set collection.}

\item{pvalues}{A list of p-values, of the same lenght as the protein list.}

\item{foldChanges}{}

\item{pianoGSC}{A piano formatted gene set collection}
}
\value{
A piano GSEA result object.
}
\description{
Runs the piano gene set enrichment analysis for a set of proteins and accompanying p-values.
}

