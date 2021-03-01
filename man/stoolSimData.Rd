\docType{data}
\name{stoolSimData}
\alias{stoolSimData}
\title{Simulated 16S rRNA gene sequencing data from stool samples}
\description{
Simulated 16S rRNA gene sequencing data from stool samples
generated with a negative binomial (NB) model (see References for further details on the simulation).
}
\usage{
stoolSimData
}
\format{
The data set contains the following objects:

    simAbund: a matrix of simulated species abundances over 8,048 species (rows) and 20 samples (columns) 
    
    simCounts: a matrix of simulated counts over 8,048 species (rows) and 20 samples (columns), resulting from the (simulated) sequencing of 'simAbund'
    
    avgAbundances: a numeric vector of average species abundances used for the simulation 
    
    phi: the coefficient of dispersion of the count data used for the simulation
    
    sdepth: a vector of sequencing depths used for the simulation 
}
\references{Measuring the diversity of the human microbiota with targeted next-generation sequencing. Briefings in Bioinformatics 19 (4), 679-692, 2018.}

\examples{
library(DiversitySeq)
data(salivaSimData)
}


