\docType{methods}
\name{divplot}
\alias{divplot}
\title{Plot diversity computed with the DiversitySeq package}
\description{
Plots the boxplots of alpha and beta diversity computed with the DiversitySeq package. 
}
\usage{
divplot(diversity, main="", ylab = "Diversity", las = 1, 
    col = NULL, points = FALSE, pointcol = "black", cexpoints = 1)
}
\arguments{
\item{diversity}{ a list of alpha or beta diversity values computed with 'adiv' or 'bdiv' functions, respectively }

\item{main}{ the title of the plot }

\item{ylab}{ the label of the y-axis }

\item{las}{ numeric in [ 0, 1, 2, 3 ] indicating the style of axis labels; for further details, check 'las' help(par) }

\item{col}{ colors to be used for the boxplots. When set to NULL, the color is set to 'white'. When set to "default", a default palette is used (up to 4 groups). Alternatively, a vector of valid color can be provided (see the examples below). }

\item{points}{ a logical value indicating whether single diversity values should be shown as scatter plot over each boxplot }

\item{pointcol}{ the color of the points in the scatter plot, when 'points = TRUE' }

\item{cexpoints}{ a numerical value giving the amount by which the points should be magnified relative to the default of 1, when 'points = TRUE'; for further details, check 'cex' in help(par) }

}

\references{Measuring the diversity of the human microbiota with targeted next-generation sequencing. Briefings in Bioinformatics 19 (4), 679-692, 2018.}

\examples{
# Load package and data set
library(DiversitySeq)
data(salivaSimData)

# Assign samples to 2 groups
group <- c(rep("group1", ncol(simCounts)/2), rep("group2", ncol(simCounts)/2))

# Compute alpha diversity with Hill numbers of order 2
adiv <- aindex(simCounts, index = "Hill", q=2, group = group)

# Plot diversity
divplot(adiv)
divplot(adiv, points = TRUE, cexpoints = 0.8)
divplot(adiv, col = "default", las = 2)
divplot(adiv, col = c("#4682B4","#FFD700"), 
        main = "Hill diversity (q=2)", 
        ylab = "Number equivalents")

}


