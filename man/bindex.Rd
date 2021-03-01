\docType{methods}
\name{bindex}
\alias{bindex}
\title{Compute beta diversity from count data}
\description{
Computes beta diversity from a matrix of species counts. 
}
\usage{
bindex(countdata, index = c("w", "c", "r", "I", 
    "e", "m", "mn", "-2","co", "cc", "-3", 
    "-3n", "rs", "sim", "z"), group = NULL)
}
\arguments{
\item{countdata}{ a matrix of species counts, with species on the rows and samples on the columns (more generally, it can be a matrix of counts computed for a set of non-overlapping classes) }

\item{index}{ the index to be used for the computation of beta diversity (see the package vignette for further details) }

\item{group}{  vector of strings indicating to which group the samples in 'countdata' belong. The length of 'group' vector must equal the number of columns of 'countdata'. When the parameter 'group' is not specified, all samples are assigned to the same group, called 'group1' }
}
\value{
List whose length equals the number of different groups. For each group, the list contains the beta diversity values computed for all pairs of samples in 'countdata' which belong to the group. 

}
\references{Measuring the diversity of the human microbiota with targeted next-generation sequencing. Briefings in Bioinformatics 19 (4), 679-692, 2018.}

\examples{
# Load package and data set
library(DiversitySeq)
data(salivaSimData)

# Assign samples to 2 groups
group <- c(rep("group1", ncol(simCounts)/2), 
rep("group2", ncol(simCounts)/2))

# Compute beta diversity with Whittaker index
bdiv <- bindex(simCounts, index = "w", group = group)

# Plot diversity
divplot(bdiv, col="default")
}


