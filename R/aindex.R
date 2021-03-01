aindex<-function(countdata, index=c("Hill", "Renyi", "BergerParker",
    "Richness", "iSimpson", "cSimpson", "Shannon",
    "Chao1", "ACE", "Jackknife1", "Jackknife2",  
    "Pielou","Tail", "EF", "IF", "RLE", "RLI"), 
    q=NULL, keep0=FALSE, scalemin=FALSE, group=NULL) {

    index<-match.arg(index, c("Hill", "Renyi","BergerParker",
        "Richness", "iSimpson", "cSimpson", "Shannon",
        "Chao1", "ACE", "Jackknife1", "Jackknife2",  
        "Pielou","Tail", "EF", "IF", "RLE", "RLI"))
    
    if (is.null(index)) {
        index<-"Richness"
    }
    
    alpha<-NA
    nc<-ncol(countdata)
    
    if (is.null(nc)) {
        
        if (!is.vector(countdata) || !is.numeric(countdata)) {
            stop("'countdata' must be a vector or matrix 
                of absolute or relative abundances\n")
        }
        
        if ((index %in% c("Hill", "Renyi", "EF", 
            "IF", "RLE", "RLI")) && is.null(q)) {
            
            stop("Please specify the order 'q'")
        }
        
        if (!is.logical(scalemin)) {
            
            stop("'scalemin' must be set to TRUE or FALSE")
        }
        
        if (!is.logical(keep0)) {
            
            stop("'keep0' must be set to TRUE or FALSE")
        }
        
        if (scalemin) {
            countdata<-countdata/min(countdata[countdata>0])
            countdata<-round(countdata)
        }
        
        # Number of present taxa
        # if keep0=TRUE #assayed=#present
        S<-qD(countdata, 0, keep0=FALSE)
        
        # Relative abundances
        K<-sum(countdata)
        pi<-countdata/K
        
        if (index == "Richness") {
            alpha<-S
            
        } else if (index == "Hill") {
            alpha<-qD(countdata, q, keep0=keep0)
            
        } else if (index == "BergerParker") {
            alpha<-1/qD(countdata, Inf, keep0=keep0)
            
        } else if (index == "Renyi") {
            alpha<-log(qD(countdata, q, keep0=keep0))
            
        } else if (index == "iSimpson") {
            alpha<-qD(countdata, 2, keep0=keep0)
            
        } else if (index == "cSimpson") {
            alpha<-1-1/(qD(countdata, 2, keep0=keep0))
            
        } else if (index == "Shannon") {
            alpha<-log(qD(countdata, 1, keep0=keep0))
            
        } else if (index == "Chao1") {
            F1<-sum(countdata==1)
            F2<-sum(countdata==2)
            if (F1==0 || F2==0) {
                alpha<-S
            } else{
                alpha<-S+(F1*(F1-1))/(2*(F2+1)) 
            }
            
        } else if (index == "Jackknife1") {
            F1<-sum(countdata==1)
            alpha<-S+F1
            
        } else if (index == "Jackknife2") {
            F1<-sum(countdata==1)
            F2<-sum(countdata==2)
            alpha<-S+2*F1-F2
            
        } else if (index == "ACE") {
            
            ACE<-estimateR(countdata)
            alpha<-ACE["S.ACE"]
            if (is.na(alpha)) alpha<-S
   
        } else if (index == "Pielou") {
            alpha<-log(qD(countdata, 1, keep0=keep0))/
                log(qD(countdata, 0, keep0=keep0))
            
        } else if (index == "Tail") {
            pi.srt<-sort(pi[pi>0], decreasing=TRUE)
            pi.srt<-pi.srt[-1]
            iv<-seq(1,length(pi.srt))
            alpha<-sqrt(sum(pi.srt*(iv^2)))
            
        } else if (index == "EF") {
            alpha<-qD(countdata, q, keep0=keep0)/
                qD(countdata, 0, keep0=keep0)
            
        } else if (index == "IF") {
            alpha<-qD(countdata, 0, keep0=keep0)/
                qD(countdata, q, keep0=keep0)
            
        } else if (index == "RLE") {
            alpha<-log(qD(countdata, q, keep0=keep0))/
                log(qD(countdata, 0, keep0=keep0))
            
        } else if (index == "RLI") {
            alpha<-1-log(qD(countdata, q, keep0=keep0))/
                log(qD(countdata, 0, keep0=keep0))
        } 
        
        # Input matrix (recursive)
        } else {
            
            if (is.null(group)) {
                group<-rep("group1", nc)
            }
            
            if (nc != length(group)) {
                stop("The length of 'group' vector must equal the number of columns of 'countdata'")     
            }
            ugroup<-unique(group)
            
            alpha<-rep(NA, nc) 
            for (cc in 1:nc) {
                alpha[cc]<-aindex(
                    as.numeric(countdata[,cc]), 
                    index=index, q=q, 
                    keep0=keep0, scalemin=scalemin)
            }
            names(alpha)<-colnames(countdata)
        }
    
    if (!is.null(nc)) {
        
        alpha.list<-vector(mode="list", length=length(ugroup))
        names(alpha.list)<-ugroup
        for (cgroup in ugroup) {
            
            alpha.list[[cgroup]]<-alpha[which(group==cgroup)]
            
        }
        alpha<-alpha.list
    } 
    
    return(alpha)
}