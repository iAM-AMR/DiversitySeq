divplot<-function(diversity, main="", ylab="Diversity", las=1, col=NULL, points=FALSE, pointcol="black", cexpoints=1) {
        
        if (!is.list(diversity)) {
                stop("'diversity' must be a list generated with 'aindex' or 'bindex' functions")
        }
        
        if (!is.logical(points)) {
                stop("'points' must be a set to TRUE or FALSE")
        }
        
        if (is.null(col)) {
                palette<-rep("white", length(diversity))
                
        } else if (length(col)==1 && col=="default") {
                
                palette<-c("#9ACD32", "#E83A5D", "#4682B4", "#FFD700")
                
                if (length(diversity) > length(palette)) {
                        warning("Too many groups for the default palette. Color set to white")
                        col<-rep("white", length(diversity))
                        
                } else {
                        
                        col<-palette[1:length(diversity)]
                }
                
        } else  {
                
                
                if (length(diversity) != length(col)) {
                        
                        warning("Wrong number of colors. Color set to white")
                        col<-rep("white", length(diversity))
                        
                } 
        }
        
        boxplot(diversity, 
                main=main, ylab=ylab, col=col, las=las)
        if (points) {
                stripchart(diversity, vertical=TRUE,
                        method = "jitter", add=TRUE, pch=20, col=pointcol, cex=cexpoints)
        }
        
}