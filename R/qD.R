qD<-function(indata, q, keep0=FALSE) {
    
    D<-NA
    
    if (!is.vector(indata) || !is.numeric(indata)) {
        stop("'indata' must be a vector of absolute 
            or relative abundances or counts")
    }
    # Relative abundances
    indata<-indata/sum(indata)
    
    tol<-1e-6
    if (abs(sum(indata)-1)>tol) {
        stop("'indata' must be a vector of absolute 
            or relative abundances or counts")
    }
    # check order
    if (length(q)!=1 || !is.numeric(q)) {
        stop("Order 'q' must be a number")
    }
    
    if (q == 0 && keep0) {  
        D<-length(indata) 
        
    } else {
        
        indata<-indata[indata>0]
        if (length(indata)<1) {
            stop("There are no present species")
        }
        
        if (q == 0) {
            D<-length(indata)
            
        } else if (q == 1) {
            p<-indata
            w<-log(indata)
            D<-exp(-(sum(p*w)))
            
        } else if (q == Inf) {
            D<-1/max(indata)
            
        } else {
            p<-indata
            w<-indata^(q-1)
            D<-(sum(p*w))^(1/(1-q))
        }
        
    }
    
    
    
    return(D)
    }
