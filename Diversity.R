# This code is to calculate Hill numbers and make figures
#Install packages below and load them

library(plyr)
library(ggplot2)
library(vegetarian)
library(reshape2)
library(patchwork)
library(directlabels)
library(viridis)

### diversity with hill number ###

com <- fernsrbcc # add here your community matrix, with species in the colunnms and sits in the rows
h <- c(0,0.5,1,1.5,2,3,5)
hilldf <- list(NA)
hillq <- matrix(NA,nrow = 7,ncol=4)
  hillq[,2]<- h

for(i in seq(nrow(com))){
  
  com2 <- com[i,]
  com3 <- decostand(com2[,which(com2>0)], method = "total", MARGIN = 1)
  
  for(q in unique(h)){
    
  if(q!=1){d <- (sum(com3^q))^(1/(1-q))
  }else {d <- exp(diversity(com3))}
    
    hillq[which(h==q),1] <-d
    hillq[,3] <- rownames(com3)
    
    e <- d/(sum(com3^0))^(1/(1-0))
    hillq[which(h==q),4] <-e
    
    #hillq[which(h==q),4] <- d/hillq[1,1]
    
    
  hilldf[[i]] <- hillq
  }
  }
  
div <- ldply(hilldf,data.frame)
colnames(div) <- c("diversity","q","ID", "eveness")

div[,1] <- as.numeric(div[,1])
  div[,1] <- round(div[,1],3)
div[,2] <- as.numeric(div[,2])
div[,3] <- as.factor(div[,3])
div[,4] <- as.numeric(div[,4])
  div[,4] <- round(div[,4],3)
  
  div2=div
  
  
  #levels(div$ID)[levels(div$ID)=="Mirador"] = "M.Sumaco"


# make somegraphics diversity 

fd <-ggplot(div,aes(x=q,y=diversity,group=ID,color=ID))+
  geom_line(size=0.8)+
  geom_dl(aes(label=ID), method=list(dl.trans(x = x + 0),"first.bumpup",cex=1))+
  theme_bw()+
  annotate(geom = "text", label= "(A)", x=5,y=40,size=7)+
  scale_y_continuous(limits = c(0, 40))+
  scale_x_continuous(limits = c(-1, 5),breaks = c(0,1,2,3,4,5))+
  theme(axis.title = element_text(size = 20))+
  theme(text = element_text(size=20))+
  theme(axis.text = element_text(color = "black"),panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  ylab("Taxonomic diversity")+
  scale_color_viridis(discrete = TRUE)+
  theme(legend.position = "none")
fd


# grafico eveness
fe <-ggplot(div,aes(x=q,y=eveness,group=ID,color=ID))+
  geom_line(size=0.8)+
  geom_dl(aes(label=ID), method=list(dl.trans(y = y + 0),"last.bumpup",cex=1))+
  theme_bw()+
  #annotate(geom = "text", label= "(b)", x=6,y=1)+
  #scale_y_continuous(limits = c(0, 40))+
  annotate(geom = "text", label= "(B)", x=6,y=1,size=7)+
  scale_x_continuous(limits = c(0, 6),breaks=c(0,1,2,3,4,5))+
  theme(axis.title = element_text(size = 20))+
  theme(text = element_text(size=20))+
  theme(axis.text = element_text(color = "black"),panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  ylab("Evenness")+
  scale_color_viridis(discrete = TRUE)+
  theme(legend.position = "none")
fe

#####################################################################
# gamma= alpha*beta

h <- c(0,0.5,1,1.5,2,3,5)
gamma <- as.data.frame(matrix(nrow = 7,ncol = 4))
colnames(gamma) <- c("a","b","g","q")

for(i in unique(h)){
  if(i==0)  {j=1}
  if(i==0.5){j=2}
  if(i==1)  {j=3}
  if(i==1.5){j=4}
  if(i==2)  {j=5}
  if(i==3)  {j=6}
  if(i==5)  {j=7}
  
  gamma[j,1] <- d(com,lev = "alpha",wts = FALSE,q=i)
  gamma[j,2] <- d(com,lev = "beta", wts = FALSE,q=i)
  gamma[j,3] <- d(com,lev = "gamma",wts = FALSE,q=i)
  gamma[4] <- h
}

gammaM <- melt(gamma,id=c("q"))
colnames(gammaM)[2] <- "index"

# grafico

pf <- ggplot(data=gammaM, aes(y=value, x=q, group=index))+
  geom_line(aes(color=index),size=1.5)+
  scale_color_grey()+
  theme_bw()+
  annotate(geom = "text", label= "(C)", x=5,y=200,size=7)+
  ylab("Diversity Partitioning")+
  scale_y_continuous(limits = c(0, 200))+
  theme(legend.position = c(0,1),legend.justification=c(-1.1,1.1),
        legend.title = element_blank(),
        legend.text = element_text(size=20))+
  theme(axis.title = element_text(size = 20))+
  theme(text = element_text(size=20))+
  theme(axis.text = element_text(color = "black"),panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  scale_color_viridis(discrete = TRUE,labels=c(expression(alpha),expression(beta),expression(gamma)))
pf


###########################################################

# save and export
jpeg("figdiv.jpg",width = 15, height = 25, units = "cm", res=300)
fd+fe+pf+ plot_layout(ncol=1,byrow = F)
dev.off()
