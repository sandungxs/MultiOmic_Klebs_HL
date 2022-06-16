
#########################################################
#           METABOLOMIC DOWNSTREAM ANALYSIS             #
#########################################################

##  Necessary packages:

#install.packages("ggplot2")
library(ggplot2, logical.return = T)
#install.packages("MetBrewer")
library(MetBrewer, logical.return = T)
#install.packages("cowplot")
library(cowplot, logical.return = T)
# install.packages("plotly")
library(plotly, logical.return = T)
# install.packages("ggrepel")
library(ggrepel, logical.return = T)
# install.packages("FactoMineR")
library(FactoMineR)
# install.packages("factoextra")
library(factoextra)


## Metabolomic data First Run:
raw.metabolomic.data.1 <- read.table(file="data/klebso_metabolomic_data_1.tsv",header=T,sep="\t",as.is=T)
head(raw.metabolomic.data.1[,3:5])
metabolomic.data.1 <- raw.metabolomic.data.1[,3:ncol(raw.metabolomic.data.1)]
weight.1 <- raw.metabolomic.data.1$Peso

# Preprocess
for(i in 1:nrow(metabolomic.data.1))
{
  print(i)
  metabolomic.data.1[i,] <- (metabolomic.data.1[i,]/weight.1[i])
}

head(metabolomic.data.1[6,])
rownames(metabolomic.data.1)<-c("LL_1","LL_2","LL_3","HL_1","HL_2","HL_3")

## Metabolomic data Second Run:
raw.metabolomic.data.2 <- read.table(file="data/klebso_metabolomic_data_2.tsv",header=T,sep="\t",as.is=T)
metabolomic.data.2 <- raw.metabolomic.data.2[,2:ncol(raw.metabolomic.data.2)]
paracetamol <- metabolomic.data.2$Paracetamol
weight.2 <- c(20.6, 17.8, 18.2, 19.6, 19.2, 17.5)

for(i in 1:nrow(metabolomic.data.2))
{
  print(i)
  metabolomic.data.2[i,] <- (metabolomic.data.2[i,]/weight.2[i])/paracetamol[i]
}

rownames(metabolomic.data.2)<-c("LL_1","LL_2","LL_3","HL_1","HL_2","HL_3")


# Getting common metabolites
length(colnames(metabolomic.data.1))
length(colnames(metabolomic.data.2))

metabolites <- intersect(colnames(metabolomic.data.1), colnames(metabolomic.data.2))
write.table(metabolites,file = "data/metabolites.txt",sep = ",")
length(metabolites)
length(setdiff(colnames(metabolomic.data.1), colnames(metabolomic.data.2))) # diferentes
length(setdiff(colnames(metabolomic.data.2), colnames(metabolomic.data.1))) # diferentes


# Data without normalization
Run_1 <- ggplot(stack(as.data.frame(t(metabolomic.data.1))), aes(x = ind, y = values, fill = ind)) +
  #stat_boxplot(geom = 'errorbar') +
  geom_boxplot(outlier.shape = NA) + 
  scale_y_continuous(limits = quantile(t(metabolomic.data.1[,2]), c(0.1, 0.8))) +
  scale_fill_manual(values=c(met.brewer("Egypt",2)[2],met.brewer("Egypt",2)[2],met.brewer("Egypt",2)[2],
                             met.brewer("Egypt",2)[1],met.brewer("Egypt",2)[1],met.brewer("Egypt",2)[1])) +
  theme(legend.position="none", axis.text.x = element_text(angle = 90), 
        plot.title = element_text(hjust = 0.5, size = 13), 
        panel.background = element_rect(fill = "white"),
        axis.title = element_text(size = 11),
        axis.text = element_text(size=9),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"),
        axis.line.x.bottom = element_line(color = 'black'),
        axis.line.y.left   = element_line(color = 'black'),
        panel.border = element_blank()
  ) +
  ggtitle("raw Data Run 1") +
  xlab("Samples") +
  ylab("Area")

Run_2 <- ggplot(stack(as.data.frame(t(metabolomic.data.2))), aes(x = ind, y = values, fill = ind)) +
  # stat_boxplot(geom = 'errorbar') +
  geom_boxplot(outlier.shape = NA) + 
  scale_y_continuous(limits = quantile(metabolomic.data.2[,3], c(0.1, 0.8))) +
  scale_fill_manual(values=c(met.brewer("Egypt",2)[2],met.brewer("Egypt",2)[2],met.brewer("Egypt",2)[2],
                              met.brewer("Egypt",2)[1],met.brewer("Egypt",2)[1],met.brewer("Egypt",2)[1])) +
  theme(legend.position="none", axis.text.x = element_text(angle = 90), 
        plot.title = element_text(hjust = 0.5, size = 13), 
        panel.background = element_rect(fill = "white"),
        axis.title = element_text(size = 11),
        axis.text = element_text(size=9),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"),
        axis.line.x.bottom = element_line(color = 'black'),
        axis.line.y.left   = element_line(color = 'black'),
        panel.border = element_blank()
  ) +
  ggtitle("Raw Data Run 2") +
  xlab("Samples") +
  ylab("Area")

plot_grid(Run_1, Run_2, labels=c("Run_1", "Run_2"), ncol = 2, nrow = 1)
ggsave2("plots/boxplot_before.png",width = 14,height = 6)


## Normalization:
# Empty needed vectors
metabolites.p.value <- vector(mode = "numeric",length=length(metabolites))
metabolites.fold.change <- vector(mode = "numeric",length=length(metabolites))
names(metabolites.p.value) <- metabolites
names(metabolites.fold.change) <- metabolites

normalize.metabolites <- matrix(nrow = length(metabolites),ncol = 12)
rownames(normalize.metabolites) <- metabolites
colnames(normalize.metabolites) <- c(paste("HL",1:6,sep="_"),paste("LL",1:6,sep="_"))

setwd(dir = "~/Escritorio/Master/TFM/shiny/figures_test/Metabolitos/Primary_Metabolites/")

# normalization and barplot generation
for(i in 1:length(metabolites))
{
  # For each metabolite
  current.metabolite <- metabolites[i]
  
  # In Run 1 for each condition
  hl.1 <- metabolomic.data.1[1:3,current.metabolite]
  ll.1 <- metabolomic.data.1[4:6,current.metabolite]
  # Standarized by the mean
  mean.ll.1 <- mean(ll.1)
  norm.hl.1 <- hl.1/mean.ll.1
  norm.ll.1 <- ll.1/mean.ll.1
  
  # In Run 2 for each condition
  hl.2 <- metabolomic.data.2[1:3,current.metabolite]
  ll.2 <- metabolomic.data.2[4:6,current.metabolite]
  # Standarized by the mean
  mean.ll.2 <- mean(ll.2)
  norm.hl.2 <- hl.2/mean.ll.2
  norm.ll.2 <- ll.2/mean.ll.2
  
  # Merging conditions
  norm.ll <- c(norm.ll.1, norm.ll.2)
  norm.hl <- c(norm.hl.1, norm.hl.2)
  
  # normalized matrix
  normalize.metabolites[i,] <- c(norm.hl,norm.ll)
  
  
  # mean and sds values
  means <- c(mean(norm.ll), mean(norm.hl))
  sds <- c(sd(norm.ll), sd(norm.hl))
  
  # If LL mean is greater than HL mean
  if(means[1] > means[2])
  {
    # Apply t Student greater than
    sig <- t.test(x = norm.ll,norm.hl,alternative = "greater")
  } else
  {
    # Apply t Student less than
    sig <- t.test(x = norm.ll,norm.hl,alternative = "less")
  }
  # Saving results
  metabolites.p.value[i] <- sig$p.value
  
  # FoldChange
  metabolites.fold.change[i] <- means[2]/means[1]
  
  # Barplot PNG
  png(filename = paste0(current.metabolite,".png"),width = 300)
  par(lwd=3)
  xpos <- barplot(means,col=c("blue","firebrick2"),
                  names.arg = c("LL","HL"),las=2,cex.names = 1.5,
                  ylim=c(0,max(means+sds)*1.3),cex.axis = 1.5,lwd=3,
                  main=current.metabolite,
                  cex.main=2)
  points(x = rep(xpos[1]+0.2,6),y=norm.ll)
  points(x = rep(xpos[2]+0.2,6),y=norm.hl)
  arrows(x0 = xpos,y0 = means + sds, x1 = xpos, y1 = means - sds,length = 0.1,
         code = 3,angle=90,lwd=2)
  
  # Significance according pvalue
  if(sig$p.value < 0.001)
  {
    points(x=xpos[2],y= (means[2] + sds[2])*1.05,pch=8)
    points(x=xpos[2]-0.1,y= (means[2] + sds[2])*1.05,pch=8)
    points(x=xpos[2]+0.1,y= (means[2] + sds[2])*1.05,pch=8)
  } else if(sig$p.value < 0.01)
  {
    points(x=xpos[2]-0.05,y= (means[2] + sds[2])*1.05,pch=8)
    points(x=xpos[2]+0.05,y= (means[2] + sds[2])*1.05,pch=8)
  } else if(sig$p.value < 0.05)
  {
    points(x=xpos[2],y= (means[2] + sds[2])*1.05,pch=8)
  }
  
  dev.off()
}


# Normality test : H0 normal, H1 no normal
pval <- c()
for (i in 1:nrow(normalize.metabolites))
{
  res <-shapiro.test(normalize.metabolites[i,])
  pval<-c(pval,res$p.value)
}
length(pval)
norm <-pval[which(pval > 0.05)]
length(norm) # de 55, 47 follow a normal distribution


# saving normalization matrixs
head(normalize.metabolites)
write.csv(normalize.metabolites,"data/normalize.metabolites.csv")
write.table(normalize.metabolites,file = "data/normalize_metabolites",sep = "\t",quote = F)
normalize.metabolites.1<-cbind(normalize.metabolites[,7:9],normalize.metabolites[,1:3]) 
head(normalize.metabolites.1)
write.table(normalize.metabolites.1,file = "data/normalize_metabolites_Run1",sep = "\t",quote = F)
normalize.metabolites.2<-cbind(normalize.metabolites[,10:12],normalize.metabolites[,4:6]) 
head(normalize.metabolites.2)
write.table(normalize.metabolites.2,file = "data/normalize_metabolites_Run2",sep = "\t",quote = F)



# Boxplot Normalized Data Run 1
Normalized.Run1 <- ggplot(stack(as.data.frame(normalize.metabolites.1)), aes(x = ind, y = values, fill = ind)) +
  #stat_boxplot(geom = 'errorbar') +
  geom_boxplot(outlier.shape = NA) + 
  scale_y_continuous(limits = quantile(normalize.metabolites.1[,2], c(0.1, 0.8))) +
  scale_fill_manual(values=c(met.brewer("Egypt",2)[2],met.brewer("Egypt",2)[2],met.brewer("Egypt",2)[2],
                             met.brewer("Egypt",2)[1],met.brewer("Egypt",2)[1],met.brewer("Egypt",2)[1])) +
  theme(legend.position="none", axis.text.x = element_text(angle = 90), 
        plot.title = element_text(hjust = 0.5, size = 13), 
        panel.background = element_rect(fill = "white"),
        axis.title = element_text(size = 11),
        axis.text = element_text(size=9),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"),
        axis.line.x.bottom = element_line(color = 'black'),
        axis.line.y.left   = element_line(color = 'black'),
        panel.border = element_blank()
  ) +
  ggtitle("Normalized Data Run 1") +
  xlab("Samples") +
  ylab("Area")

# Boxplot datos normalizados Run 2
Normalized.Run2 <- ggplot(stack(as.data.frame(normalize.metabolites.2)), aes(x = ind, y = values, fill = ind)) +
  #stat_boxplot(geom = 'errorbar') +
  geom_boxplot(outlier.shape = NA) + 
  scale_y_continuous(limits = quantile(normalize.metabolites.2[,3], c(0.1, 0.8))) +
  scale_fill_manual(values=c(met.brewer("Egypt",2)[2],met.brewer("Egypt",2)[2],met.brewer("Egypt",2)[2],
                             met.brewer("Egypt",2)[1],met.brewer("Egypt",2)[1],met.brewer("Egypt",2)[1])) +
  theme(legend.position="none", axis.text.x = element_text(angle = 90), 
        plot.title = element_text(hjust = 0.5, size = 13), 
        panel.background = element_rect(fill = "white"),
        axis.title = element_text(size = 11),
        axis.text = element_text(size=9),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"),
        axis.line.x.bottom = element_line(color = 'black'),
        axis.line.y.left   = element_line(color = 'black'),
        panel.border = element_blank()
  ) +
  ggtitle("Normalized Data Run 2") +
  xlab("Samples") +
  ylab("Area")

# Comparation of pre y post normalized data
plot_grid(Run_1, Normalized.Run1, labels=NULL, ncol = 2, nrow = 1)
plot_grid(Run_2, Normalized.Run2, labels=NULL, ncol = 2, nrow = 1)


# Interactive plot
boxplot_Run1<-ggplotly(Run_1) 
boxplot_Run2<-ggplotly(Run_2) 
boxplot_normRun1 <-ggplotly(Normalized.Run1) 
boxplot_normRun2 <-ggplotly(Normalized.Run2) 

fig <- subplot(boxplot_Run1,boxplot_normRun1) %>% 
  layout(title = 'Datos Pre y Post-Normalizados') 
fig

htmlwidgets::saveWidget(as_widget(fig), "plots/boxplot_Run1.html")

fig2 <- subplot(boxplot_Run2,boxplot_normRun2) %>% 
  layout(title = 'Datos Pre y Post-Normalizados') 
fig2

htmlwidgets::saveWidget(as_widget(fig2), "plots/boxplot_Run2.html")


# Activated and Repressed primary metabolites
activated.metabolites <- metabolites[metabolites.p.value < 0.05 & metabolites.fold.change > 1]
repressed.metabolites <- metabolites[metabolites.p.value < 0.05 & metabolites.fold.change < 1]
head(metabolites.p.value)
head(metabolites.fold.change)
head(metabolites) # 55


# Volcano plot for primary metabolite
met.results<-as.data.frame(cbind(metabolites.fold.change,metabolites.p.value))
head(met.results)
met.results[,"met_type"] <- "ns"
met.results[,"met_name"] <- rownames(met.results)

met.results[which(met.results$metabolites.fold.change > 1 & met.results$metabolites.p.value < 0.05),"met_type"] <- "activado"
met.results[which(met.results$metabolites.fold.change < 1 & met.results$metabolites.p.value < 0.05),"met_type"] <- "reprimido"

write.table(met.results,file = "Primary_Metabolites/met_results.tsv",sep = "\t") # sin carotenoides ni fitohormonas

volcol <- c(met.brewer("Egypt",3)[1], met.brewer("Egypt",3)[2],"grey33")
names(volcol) <- c("activado","reprimido","ns")


volcano_met <- ggplot(as.data.frame(met.results), aes(x=log2(metabolites.fold.change), y=-log10(metabolites.p.value),
                                                   color=met_type,
                                                   text= paste0("</br> Metabolilte: " ,met_name))) +
  geom_point(size = 1.2) +
  scale_colour_manual(values = volcol) + 
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"),
        axis.line.x.bottom = element_line(color = 'black'),
        axis.line.y.left   = element_line(color = 'black'),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 15)) + 
  ggtitle("Volcano plot from Metabolites")+ 
  labs(x = "logFC",y="adj.P.Value") +
  xlim(-2, 2) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "grey33",size=0.3)

volcano_met
ggsave2("Volcano_primary_met.png",width = 10,height = 8)

# Interactive Volcano
interactive_volcano <- ggplotly(volcano_met,x = ~logFC, y = ~adj.P.Value,
                                tooltip = "text",
                                width = 400, height = 400)
interactive_volcano
htmlwidgets::saveWidget(as_widget(interactive_volcano), "plots/interactive_volcano.html")






# PCA y clustering

colnames(normalize.metabolites)

pca.metabolites <- data.frame(colnames(normalize.metabolites[,-1]),t(normalize.metabolites[,-1]))
colnames(pca.metabolites)[1] <- "Sample"
head(pca.metabolites)

res.pca.metabolites <- PCA(pca.metabolites, graph = FALSE,scale.unit = T,quali.sup = 1)
res.hcpc.metabolites <- HCPC(res.pca.metabolites, graph=FALSE)

# HCPC

res.dist <- get_dist(pca.metabolites[,-1], stand = F, method = "spearman")
dendrograme<-hclust(res.dist)
setwd(dir = "~/Escritorio/Master/TFM/shiny/figures_test/Metabolitos/plots/")
png(filename = paste0("Dendograme.png"),width = 300)
fviz_dend(dendrograme,k=2,
          cex = 1,                       # Label size
          palette = "jco",               # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          rect_border = "jco",           # Rectangle color
          type="rectangle",
          labels_track_height = 10,      # Augment the room for labels
          ylim=c(-0.5,2.1)
          
          
)
dev.off()

# reset WD to source file location

# PCA
setwd(dir = "~/Escritorio/Master/TFM/shiny/figures_test/Metabolitos/plots/")
png(filename = paste0("PCA.png"),width = 300)
fviz_pca_ind(res.pca.metabolites, col.ind = c(rep("HL",5), rep("LL",6)), 
             pointsize=2, pointshape=21,fill="black",
             repel = TRUE, 
             addEllipses = TRUE,ellipse.type = "confidence",
             legend.title="Conditions",
             title="",
             show_legend=TRUE,show_guide=TRUE) 
dev.off() # No batch effect


# reset WD to source file location


########################
#     Carotenoids     # 
########################

# Carotenoids data
carotenoids <- read.table(file="data/carotenoids_klebsormidium.tsv",header=T,sep="\t",as.is=T)
head(carotenoids)

carotenoids.names <- carotenoids$carotenoid
carotenoids <- carotenoids[,2:ncol(carotenoids)]
rownames(carotenoids) <- carotenoids.names
head(carotenoids)
ncol(carotenoids)

# Empty needed vectors
carotenoids.p.value <- vector(mode = "numeric",length=length(carotenoids.names))
carotenoids.fold.change <- vector(mode = "numeric",length=length(carotenoids.names))
names(carotenoids.p.value) <- carotenoids.names
names(carotenoids.fold.change) <- carotenoids.names

normalize.carotenoids <- matrix(nrow = length(carotenoids.names),ncol = 7)
rownames(normalize.carotenoids) <- carotenoids.names
colnames(normalize.carotenoids) <- c(paste("HL",1:4,sep="_"),paste("LL",1:3,sep="_"))


# Normalization
setwd(dir = "~/Escritorio/Master/TFM/shiny/figures_test/Metabolitos/Carotenoids/")
for(i in 1:length(carotenoids.names))
{
  
  current.carotenoid <- carotenoids.names[i]
  
  hl.1 <- unlist(carotenoids[current.carotenoid,1:4])
  ll.1 <- unlist(carotenoids[current.carotenoid,5:7])
  
  mean.ll.1 <- mean(ll.1)
  
  norm.hl.1 <- hl.1/mean.ll.1
  norm.ll.1 <- ll.1/mean.ll.1
  
  norm.ll <- c(norm.ll.1)
  norm.hl <- c(norm.hl.1)
  normalize.carotenoids[current.carotenoid,] <- c(norm.hl,norm.ll)
  
  means <- c(mean(norm.ll), mean(norm.hl))
  sds <- c(sd(norm.ll), sd(norm.hl))
  
  if(means[1] > means[2])
  {
    sig <- t.test(x = norm.ll,norm.hl,alternative = "greater")
  } else
  {
    sig <- t.test(x = norm.ll,norm.hl,alternative = "less")
  }
  
  carotenoids.p.value[i] <- sig$p.value
  carotenoids.fold.change[i] <- means[2]/means[1]
  
  png(filename = paste0(current.carotenoid,".png"),width = 300)
  par(lwd=3)
  xpos <- barplot(means,col=c("blue","firebrick2"),
                  names.arg = c("LL","HL"),las=2,cex.names = 1.5,
                  ylim=c(0,max(means+sds)*1.2),cex.axis = 1.5,lwd=3,
                  main=current.carotenoid,
                  cex.main=2)
  points(x = rep(xpos[1]+0.2,3),y=norm.ll)
  points(x = rep(xpos[2]+0.2,4),y=norm.hl)
  arrows(x0 = xpos,y0 = means + sds, x1 = xpos, y1 = means - sds,length = 0.1,
         code = 3,angle=90,lwd=2)
  
  if(sig$p.value < 0.001)
  {
    points(x=xpos[2],y= (means[2] + sds[2])*1.05,pch=8)
    points(x=xpos[2]-0.1,y= (means[2] + sds[2])*1.05,pch=8)
    points(x=xpos[2]+0.1,y= (means[2] + sds[2])*1.05,pch=8)
  } else if(sig$p.value < 0.01)
  {
    points(x=xpos[2]-0.05,y= (means[2] + sds[2])*1.05,pch=8)
    points(x=xpos[2]+0.05,y= (means[2] + sds[2])*1.05,pch=8)
  } else if(sig$p.value < 0.05)
  {
    points(x=xpos[2],y= (means[2] + sds[2])*1.2,pch=8)
  }
  
  dev.off()
}

# reset WD to source file location


# Normality test: H0 normal, H1 no normal
pval <- c()
for (i in 1:nrow(normalize.carotenoids))
{
  res <-shapiro.test(normalize.carotenoids[i,])
  pval<-c(pval,res$p.value)
}
length(pval)
norm <-pval[which(pval > 0.05)]
length(norm) #6/9 follow a normal distribution



# normalized matrix
head(normalize.carotenoids)
write.csv(normalize.carotenoids,"data/normalize.carotenoids.csv")

# Merging with primary metabolite results
metabolites.fold.change <- c(metabolites.fold.change,carotenoids.fold.change)
metabolites.p.value <- c(metabolites.p.value,carotenoids.p.value)
metabolites <- c(metabolites,carotenoids.names)
length(metabolites) # 64 metabolites DE

# Activated and repressed
activated.metabolites_2 <- metabolites[metabolites.p.value < 0.05 & metabolites.fold.change > 1]
repressed.metabolites_2 <- metabolites[metabolites.p.value < 0.05 & metabolites.fold.change < 1]


# Volcano Plot for new data:
met.results_2<-as.data.frame(cbind(metabolites.fold.change,metabolites.p.value))
head(met.results_2)
met.results_2[,"met_type"] <- "ns"
met.results_2[,"met_name"] <- rownames(met.results_2)

met.results_2[which(met.results_2$metabolites.fold.change > 1 & met.results_2$metabolites.p.value < 0.05),"met_type"] <- "activado"
met.results_2[which(met.results_2$metabolites.fold.change < 1 & met.results_2$metabolites.p.value < 0.05),"met_type"] <- "reprimido"

volcol <- c(met.brewer("Egypt",3)[1], met.brewer("Egypt",3)[2],"grey33")
names(volcol) <- c("activado","reprimido","ns")

volcano_met_2 <- ggplot(as.data.frame(met.results_2), aes(x=log2(metabolites.fold.change), y=-log10(metabolites.p.value),
                                                      color=met_type,
                                                      text= paste0("</br> Metabolilte: " ,met_name))) +
  geom_point(size = 1.2) +
  scale_colour_manual(values = volcol) + 
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"),
        axis.line.x.bottom = element_line(color = 'black'),
        axis.line.y.left   = element_line(color = 'black'),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 15)) + 
  ggtitle("Volcano plot from Metabolites \n and carotenoids")+ 
  labs(x = "logFC",y="adj.P.Value") +
  xlim(-3, 6) +
  ylim(0,6)+
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "grey33",size=0.3)

volcano_met_2
ggsave2("plots/Volcano_carotenoids_and_primary.png",width = 10,height = 8)


# Interactive Volcano
interactive_volcano_2 <- ggplotly(volcano_met_2,x = ~logFC, y = ~adj.P.Value,
                                tooltip = "text",
                                width = 400, height = 400)
interactive_volcano_2
htmlwidgets::saveWidget(as_widget(interactive_volcano_2), "plots/interactive_volcano_carotenoids.html")


########################
#     Phytohormones    # 
########################

# Phytohormone data.
phytohormone.data <- read.table(file = "data/phytohormone.tsv",header=T)
head(phytohormone.data$Sample)
phytohormone.abundance <- phytohormone.data[,3:ncol(phytohormone.data)]
head(phytohormone.abundance)
phytohormones.name <- colnames(phytohormone.abundance)

# Standarization
for(i in 1:nrow(phytohormone.data))
{
  phytohormone.abundance[i,] <- phytohormone.abundance[i,]/phytohormone.data$Weight[i]
}

# Empty needed vectors
phytohormone.fold.change <- vector(mode = "numeric",length = 5)
names(phytohormone.fold.change) <- colnames(phytohormone.abundance)
phytohormone.p.value <- vector(mode = "numeric",length = 5)
names(phytohormone.p.value) <- colnames(phytohormone.abundance)

normalize.phytohormones <- matrix(nrow = length(phytohormones.name),ncol = 12)
rownames(normalize.phytohormones) <- phytohormones.name
colnames(normalize.phytohormones) <- c(paste("LL",1:6,sep="_"),paste("HL",1:6,sep="_"))


# Normalization and barplot generation
setwd(dir = "~/Escritorio/Master/TFM/shiny/figures_test/Metabolitos/Phytohormones/")

for(j in 1:ncol(phytohormone.abundance))
{
  current.phytohormone <- colnames(phytohormone.abundance)[j]
  
  ll.j.1.2.3 <- phytohormone.abundance[1:3,j]
  ll.j.4.5.6 <- phytohormone.abundance[4:6,j]
  hl.j.1.2.3 <- phytohormone.abundance[7:9,j]
  hl.j.4.5.6 <- phytohormone.abundance[10:12,j]
  
  norm.ll <- c(ll.j.1.2.3 / mean(ll.j.1.2.3), ll.j.4.5.6 / mean(ll.j.4.5.6))
  norm.hl <- c(hl.j.1.2.3 / ll.j.1.2.3, hl.j.4.5.6 / ll.j.4.5.6)
  normalize.phytohormones[current.phytohormone,] <- c(norm.ll,norm.hl)
  
  means <- c(mean(norm.ll), mean(norm.hl))
  sds <- c(sd(norm.ll), sd(norm.hl))
  
  if(means[1] > means[2])
  {
    sig <- t.test(x = norm.ll,norm.hl,alternative = "greater")
  } else
  {
    sig <- t.test(x = norm.ll,norm.hl,alternative = "less")
  }
  
  phytohormone.fold.change[current.phytohormone] <- mean(norm.hl)
  phytohormone.p.value[current.phytohormone] <- sig$p.value
  
  png(filename = paste0(current.phytohormone,".png"),width = 300)
  par(lwd=3)
  xpos <- barplot(means,col=c("blue","firebrick2"),
                  names.arg = c("LL","HL"),las=2,cex.names = 1.5,
                  ylim=c(0,max(means+sds)*1.2),cex.axis = 1.5,lwd=3,
                  main=current.phytohormone,
                  cex.main=2)
  points(x = rep(xpos[1]+0.2,6),y=norm.ll)
  points(x = rep(xpos[2]+0.2,6),y=norm.hl)
  arrows(x0 = xpos,y0 = means + sds, x1 = xpos, y1 = means - sds,length = 0.1,
         code = 3,angle=90,lwd=2)
  
  if(sig$p.value < 0.001)
  {
    points(x=xpos[2],y= (means[2] + sds[2])*1.2,pch=8)
    points(x=xpos[2]-0.1,y= (means[2] + sds[2])*1.2,pch=8)
    points(x=xpos[2]+0.1,y= (means[2] + sds[2])*1.2,pch=8)
  } else if(sig$p.value < 0.01)
  {
    points(x=xpos[2]-0.05,y= (means[2] + sds[2])*1.05,pch=8)
    points(x=xpos[2]+0.05,y= (means[2] + sds[2])*1.05,pch=8)
  } else if(sig$p.value < 0.05)
  {
    points(x=xpos[2],y= (means[2] + sds[2])*1.2,pch=8)
  }
  
  dev.off()
}


# Normality test: H0 normal, H1 no normal
pval <- c()
for (i in 1:nrow(normalize.phytohormones))
{
  res <-shapiro.test(normalize.phytohormones[i,])
  pval<-c(pval,res$p.value)
}
length(pval)
norm <-pval[which(pval > 0.05)]
length(norm) # 4/5 follow a normal distribution


# reset WD to source file location


# Normalized matrix
head(normalize.phytohormones)
write.csv(normalize.phytohormones,"data/normalize.phytohormones.csv")

# Adding new data
metabolites.fold.change <- c(metabolites.fold.change,phytohormone.fold.change)
metabolites.p.value <- c(metabolites.p.value,phytohormone.p.value)

length(metabolites.fold.change) # 69 metabolite DE

# Final DF
metabolites.df <- data.frame(metabolite=names(metabolites.fold.change),
                             fold.change=metabolites.fold.change,
                             p.value=metabolites.p.value)

head(metabolites.df)
metabolites.df$fold.change <-log2(metabolites.df$fold.change)

write.table(x = metabolites.df,file = "data/final_metabolite_df.tsv",row.names = F,sep = "\t")


# Activated and repressed
metabolites <- names(metabolites.fold.change)

activated.metabolites <- metabolites[metabolites.p.value < 0.05 & metabolites.fold.change > 1]
repressed.metabolites <- metabolites[metabolites.p.value < 0.05 & metabolites.fold.change < 1]
length(activated.metabolites) # 12 activated metabolite
length(repressed.metabolites) # 8 repressed metabolite

# metabolite DE/total metabolite
(length(activated.metabolites) + length(repressed.metabolites)) / length(metabolites) #29%


# Volcano Plot:
met.results_3<-as.data.frame(cbind(metabolites.fold.change,metabolites.p.value))
head(met.results_3)
length(rownames(met.results_3))
met.results_3[,"met_type"] <- "ns"
met.results_3[,"met_name"] <- rownames(met.results_3)

met.results_3[which(met.results_3$metabolites.fold.change > 1 & met.results_3$metabolites.p.value < 0.05),"met_type"] <- "activado"
met.results_3[which(met.results_3$metabolites.fold.change < 1 & met.results_3$metabolites.p.value < 0.05),"met_type"] <- "reprimido"

head(met.results_3)
write.table(met.results_3,file = "data/final_DEM_result.tsv",sep = "\t") # final

volcol <- c(met.brewer("Egypt",3)[1], met.brewer("Egypt",3)[2],"grey33")
names(volcol) <- c("activado","reprimido","ns")

volcano_met_3 <- ggplot(as.data.frame(met.results_3), aes(x=log2(metabolites.fold.change), y=-log10(metabolites.p.value),
                                                          color=met_type,
                                                          text= paste0("</br> Metabolilte: " ,met_name))) +
  geom_point(size = 4) +
  scale_colour_manual(values = volcol) + 
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"),
        axis.line.x.bottom = element_line(color = 'black'),
        axis.line.y.left   = element_line(color = 'black'),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 15),
        axis.text=element_text(size=18),
        axis.title=element_text(size=18,face="bold")) + 
  # ggtitle("Volcano plot from Metabolites \n and carotenoids and phytohormone")+ 
  labs(x = "logFC",y="adj.P.Value") +
  xlim(-5.5, 5.5) +
  ylim(0,7) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "grey33",size=0.3)+
  geom_text_repel(
    data = subset(met.results_3, met.results_3$metabolites.p.value < 0.05),
    aes(label = met_name),
    size = 7,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )

volcano_met_3
ggsave2("plots/Volcano_met.png",width = 10,height = 8)

# Interactive Volcano with all metabolites
interactive_volcano_3 <- ggplotly(volcano_met_3,x = ~logFC, y = ~adj.P.Value,
                                  tooltip = "text",
                                  width = 800, height = 600)
interactive_volcano_3
htmlwidgets::saveWidget(as_widget(interactive_volcano_3), "plots/interactive_volcano_all.html")

