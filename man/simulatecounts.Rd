\docType{methods}
\name{simulatecounts}
\alias{simulatecounts}
\title{Simulate 16S rRNA gene sequencing data}
\description{
Simulate a count matrix from 16S rRNA gene sequencing.
}
\usage{
simulatecounts(avgAbund, phi, sdepth)
}
\arguments{
\item{avgAbund}{ a numeric vector of average abundances over N species}

\item{phi}{ the coefficient of dispersion of the count data to be simulated }

\item{sdepth}{ a vector of sequencing depths for the M samples to be simulated}

}
\value{
List containing the N x M the matrix of species abundances and the N x M matrix of counts.
}
\references{Measuring the diversity of the human microbiota with targeted next-generation sequencing. Briefings in Bioinformatics 19 (4), 679-692, 2018.}

\examples{
# Load package and simulation parameters
library(DiversitySeq)
data(salivaSimData)

# Simulate a new data set
newsimdata <- simulatecounts(avgAbund = avgAbundances, phi, sdepth)
newcounts<-newsimdata$counts
newabundances<-newsimdata$abundances
}


