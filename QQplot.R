
##QQ Plots p-value
directory="/clusterdata/uqpwains/ibscratch/QC-IDATS"
setwd(directory) 
load('dmp.RObject')
p1 <- dmp$pval
p2 <- sort(t(p1))
n  <- length(p2)
k  <- c(1:n)
alpha   <- 0.05
lower   <- qbeta(alpha/2,k,n+1-k)
upper   <- qbeta((1-alpha/2),k,n+1-k)
expect  <- (k-0.5)/n
biggest <- ceiling(max(-log10(p2),-log10(expect)))
lambd = median(qchisq(1-p1,1,lower.tail = TRUE)/qchisq(0.5, df=1,lower.tail = TRUE))
shade <- function(x1, y1, x2, y2, color = col.shade) {
  n <- length(x2)
  polygon(c(x1, x2[n:1]), c(y1, y2[n:1]), border = NA, col = color)
}

#tiff(file="MND_QQ.tif", bg = "white", units = "cm", width = 16, height = 16, res = 300, compression = "lzw")
jpeg("QQ_DMPpval_nocov.jpg")
plot(-log10(expect),-log10(p2),main="PD DMR with cov adjustment - 2st try",xlim=c(0,biggest),ylim=c(0,biggest),ylab=expression(paste("Ob
served  ","-",log[10],"(P)")), xlab=expression(paste("Expected  ","-",log[10],"(P)")), type = "n")
text(-0.15, y=biggest-0.3, paste("lambda = ", signif(lambd,3),sep=""), pos = 4, font=2)
shade(-log10(expect),-log10(lower),-log10(expect),-log10(upper), color = "gray")
abline(0,1,col="white",lwd=2)
points(-log10(expect),-log10(p2), pch=20)
dev.off()




