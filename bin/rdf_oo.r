data <- read.table("goo.out", col.names=c("r", "g","n"))

r <- data$r
g <- data$g
n <- data$n

jpeg(file="rdf_oo.jpg", res=300, width=7, height=5, units='in')
par(mar=c(4.7, 5.8, 1.5,1.2))
par(mgp=c(3.7,1,0))

plot(r, g, type="lines", col="black", ylab=expression(italic(g)[O-O]), xlab=expression(italic(r)~"["~ring(A)~"]"),cex.lab=1.8,lwd=2.5, xaxs="i",yaxs="i",tck=0.02,cex.axis=1.8, las=1, xlim=c(2,5),ylim=c(0,5))
lines(r,n, type="lines", lwd=2.5, col="blue")
grid()

legend("topright", bty="n", c("RDF","Bilangan Koordinasi"), col=c("black","blue"), lwd=c(2.5,2.5))

box(which="plot")
#Turn off device driver
dev.off()

#Restore default margins
par(mar=c(5,4,4,2) + 0.1)
