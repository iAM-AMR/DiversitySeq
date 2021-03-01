mergedatasets<-function(datasets, groups) {
    
    if (!is.list(datasets)) {
        stop("'datasets' must be a list")
    }
    
    if (!is.list(groups)) {
        stop("'groups' must be a list")
    }
    
    if (length(groups)!=length(groups)) {
        stop("'datasets' and 'groups' must have the same length")
    }
    
    allspecies<-c()
    allgroups<-c()
    for (i in 1:length(datasets)) {
        
        cdata<-datasets[[i]]
        cgroup<-groups[[i]]
        
        allgroups<-c(allgroups,cgroup)
        
        if (length(cgroup)!=ncol(cdata)) {
            stop("the number of groups must be equal to the number of columns in each data set")
        }
        
        allspecies<-unique(c(allspecies, rownames(unlist(cdata))))
        
    }
    
    for (i in 1:length(datasets)) {
     
        cdata<-datasets[[i]]
        speciesnames<-rownames(cdata)
        
        if (is.null(speciesnames)) {
            stop("The row names of each data set must be not null")
        }
        
        cdata<-cdata[match(allspecies, speciesnames),]
        
        if (i==1) {
            mergedData<-cdata
            
        } else {
            mergedData<-cbind(mergedData, cdata)
            
        }
        # Set 0 counts for absent species
        mergedData[is.na(mergedData)]<-0
        
        rownames(mergedData)<-allspecies
    }
    
    return(list(data=mergedData, group=allgroups))
}