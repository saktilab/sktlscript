# Read data from the file
dataset1 <- read.table("a", header=T)
#dataset2 <- read.table("fes_data/amp/fes103.dat", header=T)
# Define colors to be used for each column
# Print pdf file
jpeg(file="pmf.jpg", res=300, width=7, height=5, units='in')
par(mar=c(4.7, 5.8, 1.5,1.2))
par(mgp=c(3.7,1,0))

x <- dataset1$x
y <- dataset1$dF*2625.5 - min(dataset1$dF*2625.5) 

#x1 <- dataset2$x
#y1 <- dataset2$dF*2625.5 - min(dataset2$dF*2625.5)

#Plot the first radial distribution function
plot(x, y, type="lines", col="black", xlab=expression(italic(n)[O-Al]), ylab=expression("Free energy"~"[kJ/mol]"),cex.lab=1.8,lwd=3.5,xlim=c(0,1), ylim=c(0,200), xaxs="i",yaxs="i",tck=0.02,cex.axis=1.8, las=1)
#lines(x1, y1, type="lines",col=2,lwd=3.5)
grid()
#lines(dataset2$r, dataset2$dE, col="red",lwd=2.5)
#Make x axis tick without labels
#axis(1,lab=F)

#Plot Other Radial Distribution Functions
#Create box around plot
box(which="plot")

#Create Legend in the top left corner that is slightly smaller and has no border
#legend("topright", c(expression("n C"~O[2]~" : n Amine = 1:2"),expression("n C"~O[2]~" : n Amine = 1:1")), lty=c(1,1,1,1,1,2,2,2,2,2), lwd=c(2.5,2.5,2.5), col=c("red", "black"), bty='n',ncol=1);

#Turn off device driver
dev.off()

#Restore default margins
par(mar=c(5,4,4,2) + 0.1)


