data <- read.table("gro.dat", col.names=c("t", "gro","index"),skip=1)

t <- data$t
gro <- data$gro

jpeg(file="gro.jpg", res=300, width=7, height=5, units='in')
par(mar=c(4.7, 5.8, 1.5,1.2))
par(mgp=c(3.7,1,0))

plot(t, gro, type="lines", col="black", ylab=expression(italic(h(t))), xlab=expression(italic(t)~"[ps]"),cex.lab=1.8,lwd=2.5, xaxs="i",yaxs="i",tck=0.02,cex.axis=1.8, las=1)

grid()

box(which="plot")
#Turn off device driver
dev.off()

#Restore default margins
par(mar=c(5,4,4,2) + 0.1)
