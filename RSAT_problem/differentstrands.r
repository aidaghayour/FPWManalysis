cebpb <- read.delim("~/Documents/R/fpwm-Thesis-new/RSAT problem/Reverse complemetn/CEBPB JASPAR format (reverse complement).ft", header=FALSE, comment.char="#")
cebpb <- data.frame(cbind(cebpb[,4],cebpb[,5], cebpb[,6],cebpb[,9]))
cebpb[,4]<-as.numeric(as.character(cebpb[,4]))
cebpb[,3]<-as.numeric(as.character(cebpb[,3]))
cebpb[,2]<-as.numeric(as.character(cebpb[,2]))
cebpb[,4]<-(-log10(cebpb[,4]))
center_pos.cebpb <- rowMeans(cebpb[,2:3])
cebpb_box <- data.frame(pos = center_pos.cebpb, pval = cebpb[,4])
binMap.cebpb <- cut( cebpb_box$pos, breaks = seq(-200,0, by = 5), labels = seq(-200,-1, by = 5 ))
boxplot(cebpb_box$pval~binMap.cebpb,ylab="-log10 Pval",xaxt="n",main="Comparision between CEBPB matrix from JASPAr and its reverse complement on RSAT",ylim=c(3,8),col = rgb(red = 1, green = 0, blue = 0, alpha = 0.3), lty=3, pch=3)


cebpb <- read.delim("~/Documents/R/fpwm-Thesis-new/RSAT problem/Reverse complemetn/CEBPB JASPAR format (Positive strand).ft", header=FALSE, comment.char="#")
cebpb <- data.frame(cbind(cebpb[,4],cebpb[,5], cebpb[,6],cebpb[,9]))
cebpb[,4]<-as.numeric(as.character(cebpb[,4]))
cebpb[,3]<-as.numeric(as.character(cebpb[,3]))
cebpb[,2]<-as.numeric(as.character(cebpb[,2]))
cebpb[,4]<-(-log10(cebpb[,4]))
center_pos.cebpb <- rowMeans(cebpb[,2:3])
cebpb_box <- data.frame(pos = center_pos.cebpb, pval = cebpb[,4])
binMap.cebpb <- cut( cebpb_box$pos, breaks = seq(-200,0, by = 5), labels = seq(-200,-1, by = 5 ))
boxplot(cebpb_box$pval~binMap.cebpb,ylab="-log10 Pval",xaxt="n",ylim=c(3,8),col = rgb(red = 0, green = 0, blue = 1, alpha = 0.3),add = TRUE, lty=1)

legend(0,7, c("Default","reverse compliment","+ : The compliment", "O : Default", ".... : The compliment", "____ : Defualt"),lty=c(1,1), lwd=c(2.5,2.5),col=c("blue","red","black","yellow","green","orange"),density = 20,cex = 0.75)


############################# Frequency

cebpb <- read.delim("~/Documents/R/fpwm-Thesis-new/RSAT problem/Reverse complemetn/CEBPB JASPAR format (reverse complement).ft", header=FALSE, comment.char="#")
reverse <- data.frame(cbind(cebpb[,4],cebpb[,5], cebpb[,6],cebpb[,9]))


cebpb <- read.delim("~/Documents/R/fpwm-Thesis-new/RSAT problem/Reverse complemetn/CEBPB JASPAR format (Positive strand).ft", header=FALSE, comment.char="#")
positive <- data.frame(cbind(cebpb[,4],cebpb[,5], cebpb[,6],cebpb[,9]))

ggplot() +geom_histogram(data = reverse,aes(X3),binwidth = 5,alpha=.2, fill="grey") +
    geom_freqpoly(data = reverse,aes(X3),binwidth = 5,col="red")+ geom_freqpoly(data = positive,aes(X3),binwidth = 5,col="blue") +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                                                        panel.background = element_blank(), axis.line = element_line(colour = "black"))+geom_density(alpha=.2, fill="#FF6666")+labs(x = "Sequence Start Point")+labs(title = "Comparison between Matrix (blue) and its reverse compliment (red)")+ geom_line(size = 2)


