bindex<-function(countdata, index=c("w", "c", "r", "I", "e", "m",
    "mn", "-2","co", "cc", "-3", "-3n", "rs", "sim", "z"), group=NULL)
{
    index<-match.arg(index, c("w", "c", "r", "I", "e", "m",
        "mn", "-2","co", "cc", "-3", "-3n", "rs", "sim", "z"))
    
    if (is.null(index)) {
        index<-"w"
    }
    
    nc<-ncol(countdata)
    if (is.null(nc) || nc<2 || !is.numeric(countdata)) {
        stop("countdata must be a matrix 
                of absolute or relative abundances\n")
    }
     
    if (is.null(group)) group<-rep("group1", nc)
    if (nc != length(group)) {
        stop("pleasy specify group for all columns of countdata\n")     
    }
    ugroup<-unique(group)
    
    beta<-vector(mode="list", length=length(ugroup))
    names(beta)<-ugroup
    for (cgroup in ugroup) {
        cdata<-countdata[,which(group==cgroup)]
        cbeta<-c()
        for (i in 1:(ncol(cdata)-1)) {
            for (j in (i+1):ncol(cdata)) {
                x<-cdata[,i]
                y<-cdata[,j]
                cbeta<-c(cbeta, bindexXY(x, y, index))
                names(cbeta)[length(cbeta)]<-paste(colnames(cdata[,c(i,j)]), collapse="_vs_")
                    
            }
        }
        beta[[cgroup]]<-cbeta
    }
     
    return(beta)
  
}