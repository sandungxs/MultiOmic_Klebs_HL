##############################
#          PAM DATA          #
##############################

# Chlorophyll fluorescence measured by Pulse-Amplitude-Modulation (PAM) Fluorometry
# for each condition in 

# Replicate 1
hl.1 <- read.table(file="HL_1.tsv",header=F)
colnames(hl.1) <- c("Time","Fluorescence")
head(hl.1)
ll.1 <- read.table(file="LL_1.tsv",header=F)
colnames(ll.1) <- c("Time","Fluorescence")
head(ll.1)

# Replicate 2
hl.2 <- read.table(file="HL_2.tsv",header=F)
colnames(hl.2) <- c("Time","Fluorescence")
head(hl.2)
ll.2 <- read.table(file="LL_2_2.csv",header=F)
colnames(ll.2) <- c("Time","Fluorescence")
head(ll.2)


# High light rep 1 representation:
plot(hl.1$Time/60,hl.1$Fluorescence,type="l",lwd=3,col="red",xlim=c(0,420/60),
     axes=F,xlab="Time (min)",ylab="Chlorophyll Fluorescence",
     cex.lab=1.5,main="High Light",cex.main=2.5)

axis(side = 1,lwd=3)
axis(side = 2,labels = c("",""),at = c(0,0.4),lwd=3)
lines(x = c(1,1),y=c(-1,1),lty=2,lwd=3)
lines(x = c(5,5),y=c(-1,1),lty=2,lwd=3)


# Low light rep 1 representation:
plot(ll.1$Time/60,ll.1$Fluorescence,type="l",lwd=3,col="blue",xlim=c(0,420/60),
     axes=F,xlab="Time (min)",ylab="Chlorophyll Fluorescence",
     cex.lab=1.5,main="Low Light",cex.main=2.5)
axis(side = 1,lwd=3)
axis(side = 2,labels = c("",""),at = c(0,0.7),lwd=3)
lines(x = c(1,1),y=c(-1,1),lty=2,lwd=3)
lines(x = c(5,5),y=c(-1,1),lty=2,lwd=3)



# High light rep 2 representation:
png(filename = "pam_high_light.png")
plot(hl.2$Time/60,hl.2$Fluorescence,type="l",lwd=3,col="red",xlim=c(0,420/60),
     axes=F,xlab="Time (min)",ylab="Chlorophyll Fluorescence",
     cex.lab=1.5,main="High Light",cex.main=2.5)

axis(side = 1,lwd=3)
axis(side = 2,labels = c("",""),at = c(0,0.4),lwd=3)
lines(x = c(1,1),y=c(-1,1),lty=2,lwd=3)
lines(x = c(5,5),y=c(-1,1),lty=2,lwd=3)
dev.off()


# Low light rep 2 representation:
png(filename = "pam_low_light.png")
plot(ll.2$Time/60,ll.2$Fluorescence,type="l",lwd=3,col="blue",xlim=c(0,420/60),
     axes=F,xlab="Time (min)",ylab="Chlorophyll Fluorescence",
     cex.lab=1.5,main="Low Light",cex.main=2.5)
axis(side = 1,lwd=3)
axis(side = 2,labels = c("",""),at = c(0,0.7),lwd=3)
lines(x = c(1,1),y=c(-1,1),lty=2,lwd=3)
lines(x = c(5,5),y=c(-1,1),lty=2,lwd=3)
dev.off()


