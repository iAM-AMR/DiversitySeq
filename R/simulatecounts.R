simulatecounts<-function(avgAbund, phi, sdepth) {
    
    # Check inputs
    if (!is.vector(avgAbund) || !is.numeric(avgAbund) || !is.null(dim(avgAbund))) {
        
        stop("'avgAbund' must be a numeric vector of species average abundances")
    }
    if (!is.vector(sdepth) || !is.numeric(sdepth) || !is.null(dim(sdepth))) {
        
        stop("'sdepth' must be a numeric vector of sample sequencing depths")
    }
    if (!is.numeric(phi) || length(phi)!=1) {
        
        stop("'phi' must be a single numeric value representing the coefficient of dispersion")
    }
    
    # Number of samples
    nsamples<-length(sdepth)
    # Number of species
    nspecies<-length(avgAbund)
    message("Simulating a count data set with ", nsamples, 
            " samples and ", nspecies, " species...")    
    
    # Relative abundances
    avg_rel_abund<-avgAbund/sum(avgAbund)
    
    # Simulate biological variability due to sampling (Gamma)
    gamma.mean<-avg_rel_abund
    gamma.var<-phi*avg_rel_abund^2
    gamma.shape<-1/phi
    gamma.scale<-phi*gamma.mean
    SimAbundances<-matrix(0,nrow=nspecies,ncol=nsamples)
    for (i in 1:nspecies)
    {
        SimAbundances[i,]<-rgamma(nsamples,
            shape=gamma.shape,
            scale=gamma.scale[i])
        
    }
    
    # No negative abundances
    SimAbundances[SimAbundances<0]<-0
    
    # Compute absolute abundances
    SimRelAbundances<-SimAbundances
    SimAbundances<-round(SimAbundances*sum(avgAbund))
    
    # Simulate technical variability due to sequencing (Poisson)
    SimCounts<-matrix(0,nrow=nspecies,ncol=nsamples)
    for (i in 1:nrow(SimRelAbundances)) {
        
        for (j in 1:ncol(SimRelAbundances)) {
            
            SimCounts[i,j]<-round(rpois(1,SimRelAbundances[i,j]*sdepth[j])) 
        }
    }
    
    # Add species names
    rownames(SimCounts)<-rownames(SimAbundances)<-names(avgAbund)
    colnames(SimCounts)<-colnames(SimAbundances)<-paste("Sample", seq(1:ncol(SimCounts)), sep="_")
        
    return(list(abundances=SimAbundances,
        counts=SimCounts))
}