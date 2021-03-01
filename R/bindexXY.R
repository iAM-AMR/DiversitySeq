bindexXY<-function(x, y, index=c("w", "c", "r", "I", "e", "m",
    "mn", "-2","co", "cc", "-3", "-3n", "rs", "sim", "z"))
{
    
    index<-match.arg(index, c("w", "c", "r", "I", "e", "m",
        "mn", "-2","co", "cc", "-3", "-3n", "rs", "sim", "z"))
    
    if (is.null(index)) {
        index<-"w"
    }
    
    beta<-NA
    
    if (!is.vector(x) || !is.numeric(x)) {
        stop("'x' must be a vector or matrix of absolute
            or relative abundances")
    }
    if (!is.vector(y) || !is.numeric(y)) {
        stop("'y' must be a vector or matrix of absolute
            or relative abundances")
    }
    
    # Shared/unique species
    xp<-which(x>0)
    yp<-which(y>0)
    a<-length(intersect(xp,yp))
    b<-length(setdiff(xp,yp))
    c<-length(setdiff(yp,xp))
    
    if (index=="w") { 
        beta<-(b+c)/(2*a+b+c)
        
    } else if (index=="c") {
        beta<-(b+c)/2
        
    } else if (index=="r") {
        beta<-2*b*c/((a+b+c)^2-2*b*c)
        
    } else if (index=="I") {
        beta<- log(2*a+b+c) - 2*a*log(2)/(2*a+b+c) - 
            ((a+b)*log(a+b) + (a+c)*log(a+c)) / (2*a+b+c)
        
    } else if (index=="e") {
        beta<-exp(log(2*a+b+c) - 2*a*log(2)/(2*a+b+c) - 
                ((a+b)*log(a+b) + (a+c)*log(a+c)) / (2*a+b+c))-1
        
    } else if (index=="m") {
        beta<-(2*a+b+c)*(b+c)/(a+b+c)
        
    } else if (index=="mn") {
        beta<-(2*a+b+c)*(b+c)/(a+b+c)^2
        
    } else if (index=="-2") {
        beta<-pmin(b,c)/(pmax(b,c)+a)
        
    } else if (index=="co") {
        beta<-(a*c+a*b+2*b*c)/(2*(a+b)*(a+c))
        
    } else if (index=="cc") {
        beta<-(b+c)/(a+b+c)
        
    } else if (index=="-3") {
        beta<-pmin(b,c)/(a+b+c)
        
    } else if (index=="-3n") {
        beta<-2*pmin(b,c)/(a+b+c)
        
    } else if (index=="rs") {
        beta<-2*(b*c+1)/((a+b+c)^2-(a+b+c))
        
    } else if (index=="sim") {
        beta<-pmin(b,c)/(pmin(b,c)+a)
        
    } else if (index=="z") {
        beta<-(log(2)-log(2*a+b+c)+log(a+b+c))/log(2)
        
    }
    
    return(beta)
    
}
