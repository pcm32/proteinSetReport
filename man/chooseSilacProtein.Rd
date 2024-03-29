% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/proteinSetReport.R
\name{chooseSilacProtein}
\alias{chooseSilacProtein}
\title{chooseSilacProtein}
\usage{
chooseSilacProtein(dataTable, colForProtsAccessions = "Protein_IDs",
  colForPeptideCounts = "Peptide_counts_unique",
  nameForNewCol = "Blessed_Protein_ID", separator = ";")
}
\arguments{
\item{dataTable}{The SILAC data set.}

\item{colForProtsAccessions}{The column name where the protein accessions are found. Defaults to "Protein_IDs".}

\item{colForPeptideCounts}{The column name where peptide counts per protein are found. Defaults to "Peptide_counts_unique".}

\item{nameForNewCol}{The name for the new column that will contain the chosen accession. Defaults to "Blessed_Protein_ID".}

\item{separator}{The token separating proteins and counts in their respective column. Defaults to ";"}
}
\value{
Same provided dataTable with the additional column.
}
\description{
Given the multiple proteins available per line on a SILAC result, choose among those with the higher amount o
f peptides,
the one that is best annotated (Swissprot vs TrEMBL)
}

