\name{DivesrsitySeq}
\alias{DivesristySeq}
\docType{package}
\title{
DiversitySeq: measuring diversity from count data sets
}
\description{
Next-generation sequencing, and particularly 16S ribosomal RNA (16S rRNA) gene sequencing, is a powerful technique for the identification and quantification of human-resident microbes, collectively known as the human microbiota. 

Once bacterial abundances are profiled via 16S rRNA gene sequencing and summarized in a count data set, diversity indices provide valuable mathematical tools to investigate the composition of the human microbiota. In brief, alpha diversity can be used to describe the compositional complexity of a single sample, whereas beta diversity can be used to identify taxonomical differences between samples. 

The DiversitySeq package implements in a unified framework the whole panel of diversity indices reviewed in Finotello et al. (2016), enabling the assessment of diversity from count data sets. DiversitySeq also implements a simulator for the generation of synthetic count data sets from 16S rRNA gene sequencing. 

Besides 16S rRNA gene sequencing data, this package can be employed with other data sets with similar characteristics, such as 5S rRNA gene sequencing, environmental metagenomics or, more generally, any kind of matrix were counts are computed for different types non-overlapping classes.
}
\details{
\tabular{ll}{
Package: \tab DiversitySeq\cr
Type: \tab Package\cr
Version: \tab \Sexpr{packageVersion("DiversitySeq")} \cr
Date: \tab 2016-01-09\cr
License: \tab GPL-3\cr
}
}
\author{
Francesca Finotello, Eleonora Mastrorilli, Barbara Di Camillo

Maintainer: Francesca Finotello <francesca.finotello@dei.unipd.it>
}
\references{Measuring the diversity of the human microbiota with targeted next-generation sequencing. Briefings in Bioinformatics 19 (4), 679-692, 2018.}

\keyword{ diversity
sequencing
NGS 
microbiota
microbiome
counts
bacteria
16S ribosomal RNA (rRNA) gene
evenness
richness
alpha
beta
}
\examples{
# Load package and data
library(DiversitySeq)
data(salivaSimData)

# Compute alpha diversity
alphadiv <- aindex(simCounts, index = "Richness")
divplot(alphadiv)

# Compute beta diversity
betadiv <- bindex(simCounts, index = "r")

# Plot beta diversity
divplot(betadiv)

# Simulate new count data
newsimdata <- simulatecounts(avgAbund = avgAbundances, phi, sdepth)
}
