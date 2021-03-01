\docType{methods}
\name{aindex}
\alias{aindex}
\title{Compute alpha diversity from count data}
\description{
Computes alpha diversity from a matrix of species counts. 
}
\usage{
aindex(countdata, index = c("Hill", "Renyi", "BergerParker",
    "Richness", "iSimpson", "cSimpson", "Shannon", 
    "Chao1", "ACE", "Jackknife1", "Jackknife2",  
    "Pielou","Tail", "EF", "IF", "RLE", "RLI"), 
     q = NULL, keep0 = FALSE, scalemin = FALSE, 
     group = NULL)
}
\arguments{
\item{countdata}{ a matrix of species counts, with species on the rows and samples on the columns (more generally, it can be a matrix of counts computed for a set of non-overlapping classes) }

\item{index}{ the index to be used for the computation of alpha diversity (see the package vignette for further details) }

\item{q}{ a number indicating the order of diversity, mandatory for Hill, Renyi, EF, IF, RLE and RLI indices. }

\item{keep0}{ a logical value TRUE/FALSE indicating whether species with null counts should be considered in the computation of alpha diversity }

\item{scalemin}{ a logical value indicating whether count data shold be scaled so to have the minimum equal to 1 (useful in case of normalized data)}

\item{group}{ a vector of strings indicating to which group the samples in 'countdata' belong. The length of 'group' vector must equal the number of columns of 'countdata'. When the parameter 'group' is not specified, all samples are assigned to the same group, called 'group1' }

}
\value{
List whose length equals the number of different groups. For each group, the list contains the alpha diversity values computed for all samples in 'countdata' which belong to the group. 

}
\references{Measuring the diversity of the human microbiota with targeted next-generation sequencing. Briefings in Bioinformatics 19 (4), 679-692, 2018.}

\examples{
# Load package and data set
library(DiversitySeq)
data(salivaSimData)

# Assign samples to 2 groups
group <- c(rep("group1", ncol(simCounts)/2), 
rep("group2", ncol(simCounts)/2))

# Compute alpha diversity with Hill numbers of order 2
adiv <- aindex(simCounts, index = "Hill", q=2, group = group)

# Plot diversity
divplot(adiv, col="default")
}


