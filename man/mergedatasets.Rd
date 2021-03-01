\docType{methods}
\name{mergedatasets}
\alias{mergedatasets}
\title{Merge count data sets}
\description{
Merge count data sets.
}
\usage{
mergedatasets(datasets, groups)
}
\arguments{
\item{datasets}{ a list containing the count matrices to be merged. Each count matrices must be a matrix of species counts, with species on the rows and samples on the columns (or, more generally, a matrix of counts conputed for a set of non-overlapping classes) }

\item{groups}{ a list containing the group annotation for each data set (in the same order as in the list above). Each 'group' vector is a vector of strings indicating to which group the samples in countdata belong. The length of 'group' vector must equal the number of columns of the corresponding count matrix }
}
\value{
List containing the merged matrix and the the corresponding vector of groups. The number of rows in the output matrix is equal to the union of all species assayed in the input data sets, and the number of columns is the sum of the samples of all input data sets.
}
\references{Measuring the diversity of the human microbiota with targeted next-generation sequencing. Briefings in Bioinformatics 19 (4), 679-692, 2018.}

\examples{
# Load package and data sets
library(DiversitySeq)
data(salivaSimData)
saliva.counts <- simCounts
data(stoolSimData)
stool.counts <- simCounts

# Generate 'group' vectors for the two data sets
saliva.group <- rep("Saliva", ncol(saliva.counts))
stool.group <- rep("Stool", ncol(stool.counts))

# Merge the data
mrgData <- mergedatasets(list(stool.counts, saliva.counts),
                        list(stool.group, saliva.group))
mrg.counts <- mrgData$data
mrg.group <- mrgData$group

}

