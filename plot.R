#Read the data
data = read.table("results.txt")
s = data$V3

#Set size of the lattice and MC steps (check in code)
N = 25
steps = 300

#Create a position vector
x = c(1:N)
y = x

#Set (plotting) size of each node
nodesize = 2.3

#Add 1 for the initial conditions
steps = steps + 1

#Create a directory to save the frames
dir.create("./frames")

for(i in 1:steps) {
  #Create a new .PNG file
  png(paste("./frames/",i,".png",sep=""),height=520,width=500)
  
  #Create a plot and set desired graphics
  par(cex.lab=1.3,cex.axis=1.3,cex.main=1.6,xpd=NA)
  plot(0,0,type="n",xlim=c(0,N+1),ylim=c(0,N+1),asp=1,xlab="x",ylab="y",main=paste("Metropolis algorithm: ",N,"×",N," lattice.",sep=""))
  box(col="black")
  legend("top", cex=1, inset=c(0,-0.065), bty="n", horiz=TRUE, legend=c("spin +1    ","spin -1"), pch=c(15,15), col=c("#F2CBE0","#C6D8FF"))
  
  for(j in 1:N) {
    for(k in 1:N) {
      if(s[(i-1)*N^2+(j-1)*N+k] == 1) lines(x[j],y[k],type="p",pch=15,col="#F2CBE0",cex=nodesize)
      if(s[(i-1)*N^2+(j-1)*N+k] == -1) lines(x[j],y[k],type="p",pch=15,col="#C6D8FF",cex=nodesize)
    }
  }
  dev.off()
}