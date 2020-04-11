#TO BE SET BY THE USER: N (number of nodes PER SIDE, not total), frames to generate,
#title of graph, file to read, and color for each neuron state. Notice that number of
#frames =/= Monte Carlo steps. This script reads every MC step calculated by the program,
#but only plots as many frames as required for a better performance. For test runs, set a 
#low number of frames to check quickly the evolution of the system. 

N = 500
frames = 30
title = "Ising model: 500x500 lattice."
filename = "./results.txt"
col1 = "black"
col2 = "orange"

#FROM HERE ON, THE SCRIPT DOES NOT NEED INTERACTION. 

#Read the data
s = read.table(filename)$V1

#Compute number of Monte Carlo steps
steps = length(s)/N^2

#Set node size (for a 500x500 px .png file)
nodesize = 57/N

#Create position vectors
x = y = c(1:N)

#Create a directory to save the frames
dir.create("./frames")

for(i in 1:steps) {
  if(i%%round(steps/frames,0) == 0) {
  #Create a new .PNG file
  png(paste("./frames/",i,".png",sep=""),height=500,width=500)
  
  #Create a plot and set desired graphics
  par(cex.lab=1.3,cex.axis=1.3,cex.main=1.6,xpd=NA)
  plot(0,0,type="n",xlim=c(0,N+1),ylim=c(0,N+1),asp=1,xlab="x",ylab="y",main=title)
  box(col="black")
  legend("top", cex=1, inset=c(0,-0.065), bty="n", horiz=TRUE, legend=c("spin +1    ","spin -1"), pch=c(15,15), col=c(col1,col2))
  
  for(j in 1:N) {
    for(k in 1:N) {
      if(s[(i-1)*N^2+(j-1)*N+k] == 1) lines(x[j],y[k],type="p",pch=15,col=col1,cex=nodesize)
      if(s[(i-1)*N^2+(j-1)*N+k] == -1) lines(x[j],y[k],type="p",pch=15,col=col2,cex=nodesize)
    }
  }
  dev.off()
  }
}