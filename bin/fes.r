data <- read.table("fes.dat", col.names=c("cv", "F"))

cv <- data$cv
F <- data$F*2625.5/4.184

jpeg(file="fes.jpg", res=300, width=7, height=5, units='in')
par(mar=c(4.7, 5.8, 1.5,1.2))
par(mgp=c(3.7,1,0))

plot(cv, F, type="lines", col="black", xlab=expression(italic(n)[Ce-O]), ylab=expression("Free energy"~"[kcal/mol]"),cex.lab=1.8,lwd=2.5,xlim=c(6.5,7), xaxs="i",yaxs="i",tck=0.02,cex.axis=1.8, las=1, c(-110,-60))

grid()

box(which="plot")
#Turn off device driver
dev.off()

#Restore default margins
par(mar=c(5,4,4,2) + 0.1)
