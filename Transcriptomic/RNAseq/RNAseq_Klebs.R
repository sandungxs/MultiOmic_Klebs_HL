#########################################
##    RNASEQ DOWNSTREAM ANALYSIS        #
#########################################

##  Necessary packages:

# For the different representations of the data, functions from the ggplot2 package
# are mainly used, and one of the color palettes interpretable by people with color blindness
# offered in the MetBrewer package. 
# install.packages("ggplot2")
library(ggplot2, logical.return = T)
#install.packages("MetBrewer")
library(MetBrewer, logical.return = T)
#install.packages("cowplot")
library(cowplot, logical.return = T)
# install.packages("plotly")
library(plotly, logical.return = T)
# install.packages("tidyr")
library(tidyr, logical.return = T) 

# Transcriptomic data of our study:
gene.expression_total <- read.table(file = "gene_expression.tsv",header=T,as.is=T)
head(gene.expression_total)
gene.ids <- gene.expression_total$geneID
gene.expression <- as.matrix(gene.expression_total[,2:ncol(gene.expression_total)])
rownames(gene.expression) <- gene.ids
head(gene.expression) # data without geneID

# Although MARACAs perform a complete downstream analysis to verify this results, 
# we perform the same downstream analysis manually.


# We use a scatterplot function to perform a primary exploratoy analysis to check
# sample variance between replicates and conditions.

scatterplot_function<-function(x,y)
{
  ggplot(as.data.frame(log2(gene.expression)+1), 
         aes(x=log2(gene.expression[,x])+1, y=log2(gene.expression[,y])+1)) + 
    geom_point(size=1, color=met.brewer("Cassatt2",1) [1],size=0.4,
               aes(text= paste0("</br> X: ",round(log2(gene.expression[,x])+1,digits = 5),
                                "</br> Y: ", round(log2(gene.expression[,y])+1,digits = 5)))) +
    xlim(0,NA) +
    ylim(0,NA) +
    geom_smooth(method=lm,color="darkgreen",size=0.3) +
    theme(panel.background = element_rect(fill = "white"),
          panel.grid.major = element_line(colour = "white"),
          panel.grid.minor = element_line(colour = "white"),
          axis.line.x.bottom = element_line(color = 'black'),
          axis.line.y.left   = element_line(color = 'black'),
          panel.border = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 13)) +
    ggtitle("Scatterplot") +
    annotate(geom = "text", x = 1.5, y = max(c(log2(gene.expression[,x]+1),log2(gene.expression[,y]+1)))-2, 
             label = paste(c(round(100*cor(gene.expression[,x],
                                           gene.expression[,y]),
                                   digits = 2),
                             "%"), collapse=""),size=2.5)
}


A_1<-ggplotly(scatterplot_function(1,1),tooltip = "text")
A_2<-ggplotly(scatterplot_function(1,2),tooltip = "text")
A_3<-ggplotly(scatterplot_function(1,3),tooltip = "text")
A_4<-ggplotly(scatterplot_function(1,4),tooltip = "text")

B_1<-ggplotly(scatterplot_function(2,1),tooltip = "text")
B_2<-ggplotly(scatterplot_function(2,2),tooltip = "text")
B_3<-ggplotly(scatterplot_function(2,3),tooltip = "text")
B_4<-ggplotly(scatterplot_function(2,4),tooltip = "text")

C_1<-ggplotly(scatterplot_function(3,1),tooltip = "text")
C_2<-ggplotly(scatterplot_function(3,2),tooltip = "text")
C_3<-ggplotly(scatterplot_function(3,3),tooltip = "text")
C_4<-ggplotly(scatterplot_function(3,4),tooltip = "text")

D_1<-ggplotly(scatterplot_function(4,1),tooltip = "text")
D_2<-ggplotly(scatterplot_function(4,2),tooltip = "text")
D_3<-ggplotly(scatterplot_function(4,3),tooltip = "text")
D_4<-ggplotly(scatterplot_function(4,4),tooltip = "text")

subplot(A_1,A_2,A_3,A_4,B_1,B_2,B_3,B_4,C_1,C_2,C_3,C_4,D_1,D_2,D_3,D_4,
        nrows = 4,shareX = F,shareY = F)



## Apply upper quantile normalization:
upper.quantiles <- vector(mode="numeric",length=ncol(gene.expression))

for(i in 1:ncol(gene.expression))
{
  upper.quantiles[i] <- quantile(gene.expression[,i],probs=0.75)
}

mean.upper.quantiles <- mean(upper.quantiles)

for(i in 1:ncol(gene.expression))
{
  gene.expression[,i] <- (gene.expression[,i] / upper.quantiles[i]) * mean.upper.quantiles
}

## Log2 transformation
log.gene.expression <- log2(gene.expression+1)
head(log.gene.expression)
write.table(log.gene.expression, file = "log_gene-expression.tsv",sep = "\t",quote = F )

# To check the normalization we permorf boxplots (before and after normalization):

A <- ggplot(stack(as.data.frame(gene.expression)), aes(x = ind, y = values, fill = ind)) +
  #stat_boxplot(geom = 'errorbar') +
  geom_boxplot(outlier.shape = NA) + 
  scale_y_continuous(limits = quantile(gene.expression[,1], c(0.1, 0.8))) +
  scale_fill_manual(values=c(met.brewer("Egypt",2)[1],met.brewer("Egypt",2)[1],
                             met.brewer("Egypt",2)[2],met.brewer("Egypt",2)[2])) +
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
  ggtitle("Raw data") +
  xlab("samples") +
  ylab("FPKM")


B <- ggplot(stack(as.data.frame(log.gene.expression)), aes(x = ind, y = values, fill = ind)) +
  #stat_boxplot(geom = 'errorbar') +
  geom_boxplot(outlier.shape = NA) + 
  scale_y_continuous(limits = quantile(log.gene.expression[,1], c(0.1, 0.8))) +
  scale_fill_manual(values=c(met.brewer("Egypt",2)[1],met.brewer("Egypt",2)[1],
                             met.brewer("Egypt",2)[2],met.brewer("Egypt",2)[2])) +
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
  ggtitle("Normalized data") +
  xlab("Samples") +
  ylab("FPKM")

# plot-grid function (cowplot packages) allows to represent several ggplots in an unique plot,
# like par(mfrow).
plot_grid(A, B, labels=c("A", "B"), ncol = 2, nrow = 1)

# To make them interactive we use the function ggplotly
boxplot_A<-ggplotly(A) 
boxplot_B<-ggplotly(B) 

# To remove visible outliers
boxplot_A$x$data <- lapply(boxplot_A$x$data, FUN = function(x){
  # Adding transparency
  x$marker = list(opacity = 0)
  # Making them smaller
  #x$marker = list(size = 0.001)
  # With the cursor you can only see the maximum value of outliers, but not each outlier.
  x$hoverinfo= list("y"< quantile(gene.expression[,1],probs = 0.75))
  return(x)
})

boxplot_A

# To remove visible outliers
boxplot_B$x$data <- lapply(boxplot_B$x$data, FUN = function(x){
  # Adding transparency
  x$marker = list(opacity = 0)
  # Making them smaller
  #x$marker = list(size = 0.001)
  # With the cursor you can only see the maximum value of outliers, but not each outlier.
  x$hoverinfo= list("y"< quantile(gene.expression[,1],probs = 0.75))
  return(x)
})

boxplot_B

# Plot interactivo con los dos
fig <- subplot(boxplot_A, boxplot_B) %>% 
  layout(title = 'Data Pre y Post-Normalization') 
fig


# To study similaritys between all conditions at the same time, we performed
# a principal component analysis (PCA) in 3D:

# We get the principal components of our gene.expresion matriz:
pca_all_samples <- prcomp(x=t(log.gene.expression),center = T,scale. = F)
PCs <- as.data.frame(pca_all_samples$x)
head(PCs)

# Defining number of PC to study
x_axis <- 1
y_axis <- 2
z_axis <- 3

# Created a completed dataframe
targets <- colnames(log.gene.expression) # samples
targets <- data.frame(targets) # samples como df
targets <- separate(targets, col = targets, sep = '_', into = c('condition','replicate')) # Separamos la condicion de la replica
targets$sample <- colnames(log.gene.expression) # columna con la muestra completa
head(targets)

# Color of each condition
color_variable <- targets$condition
# variable with sample names
id <- targets$sample


# Axis name
xlab <- paste0('PC',x_axis,'\n(', round((pca_all_samples$sdev^2)[x_axis]/sum(pca_all_samples$sdev^2)*100, 2), '%)')
ylab <- paste0('PC',y_axis,'\n(', round((pca_all_samples$sdev^2)[y_axis]/sum(pca_all_samples$sdev^2)*100, 2), '%)')
zlab <- paste0('PC',z_axis,'\n(', round((pca_all_samples$sdev^2)[z_axis]/sum(pca_all_samples$sdev^2)*100, 2), '%)')

# Plotly PCA:
trace <- list()
i=1

# Color information
for (x in unique(color_variable)){
  
  trace[[i]] <- list(mode="markers", 
                     name = x,
                     type = "scatter3d",
                     x = PCs[,1][color_variable==x],
                     y = PCs[,2][color_variable==x],
                     z = PCs[,3][color_variable==x],
                     text = id[color_variable==x])
  i <- i + 1
}



# Plot style
layout <- list(
  scene = list(
    xaxis = list(
      title = xlab, 
      showline = FALSE), 
    yaxis = list(
      title = ylab, 
      showline = FALSE), 
    zaxis = list(
      title = zlab, 
      showline = FALSE)), 
  title = "PCA (3D)")

# Final PCA 3D plot:
p <- plot_ly()

for (x in trace)
  {
  p <- add_trace(p, mode=x$mode, name=x$name, type=x$type, x=x$x, y=x$y, z=x$z, text=x$text)
  }
# Rellenamos con el estilo del plot
p <- layout(p, scene=layout$scene, title=layout$title) # only observable with a browser
                                                       # (google chrome, for example)

# Save in html format.
setwd('~/Escritorio/Master/TFM/Datos_Artículo/knitens_rnaseq/RNAseq/') # salida del html
htmlwidgets::saveWidget(as_widget(p), "PCA3D.html", selfcontained = T)



## Mean expression of MARACAs result (without normalization) for our study:
gene.expression <- as.data.frame(gene.expression)
gene.expression.ll <- (gene.expression$LL_1 + gene.expression$LL_2)/2
gene.expression.hl <- (gene.expression$HL_1 + gene.expression$HL_2)/2

mean.gene.expression <- data.frame(geneID=rownames(gene.expression),
                                   LL=gene.expression.ll,
                                   HL=gene.expression.hl)
head(mean.gene.expression)



library(limma)

## Specification of the experimental design
limma.experimental.design <- model.matrix(~ -1+factor(c(1,1,2,2)))
colnames(limma.experimental.design) <- c("LL","HL")

## Linear model fit
linear.fit <- lmFit(log.gene.expression, limma.experimental.design)

## Contrast specification and computation
contrast.matrix <- makeContrasts(contrasts = "HL-LL",
                                 levels=c("LL","HL"))

contrast.linear.fit <- contrasts.fit(linear.fit, contrast.matrix)
contrast.results <- eBayes(contrast.linear.fit)

## Extract results
de.results <- topTable(contrast.results, number=nrow(gene.expression),coef=1,sort.by="logFC")

# Activated and Repressed gene determination:
fold.change <- de.results$logFC
q.values <- de.results$adj.P.Val
genes.ids <- rownames(de.results)

names(fold.change) <- genes.ids
names(q.values) <- genes.ids

activated.genes <- genes.ids[fold.change > log2(2) & q.values < 0.05]
repressed.genes <- genes.ids[fold.change < - log2(2) & q.values < 0.05]

length(activated.genes) # 677 genes
length(repressed.genes)  # 678 genes


# This values match with MARACAs result.
act <-read.table("activated_genes.txt",header = F) # 677 genes
length(act$V1)
rep <-read.table("repressed_genes.txt",header = F) # 678 genes
length(rep$V1)







# Volcano ggplot
library(ggrepel) # Adding label without overlap

# Variable for expression level
de.results[,"gene_type"] <- "ns" 
head(de.results)
de.results[,"gene_name"] <- rownames(de.results)

# Se modifica ns por activado o reprimido:
de.results[which(de.results$logFC > 1 & de.results$adj.P.Val < 0.05),"gene_type"] <- "activated"
de.results[which(de.results$logFC< -1 & de.results$adj.P.Val < 0.05),"gene_type"] <- "repressed"

genes.ids.de.results <- rownames(de.results)

# Adding annotation.
anotacion<-read.table("annotation_klebs.txt",header = T,sep = "\t")
head(anotacion)
de.results[,"annotation"]<-"hypothetical protein"
head(de.results,20)
nrow(de.results)

# Adding annotation.
for (i in 1:nrow(anotacion))
{
  coincidencia<-subset(anotacion,anotacion[,1] == de.results[i,"gene_name"])
  if(is.na(coincidencia[1,3])== T)
  {
    de.results[i,"annotation"]<-"hypothetical protein"
  } else {
    de.results[i,"annotation"]<-coincidencia[1,3]
  }
}


# Merging gene.expression with de.results
resultado_final<-matrix(ncol = 14)
colnames(resultado_final)<-c(colnames(gene.expression_total),colnames(de.results))

for (i in 1:nrow(gene.expression_total))
{
  coincidencia<-subset(gene.expression_total,gene.expression_total[,1] == de.results[i,"gene_name"])
  fila<-cbind(coincidencia, de.results[i,])
  resultado_final<-rbind(resultado_final, fila)
}

resultado_final<-resultado_final[-1,]
rownames(resultado_final)<-NULL
head(resultado_final)

write.table(resultado_final, file = "resultado_final.tsv",sep = "\t")
# head(read.table(file = "resultado_final.tsv",header=T,as.is=T))


# Volcano plot colors.
volcol <- c(met.brewer("Egypt",3)[1], met.brewer("Egypt",3)[2],"grey33")

# Assign colors by gene expression level.
names(volcol) <- c("activated","repressed","ns")

volcano_RNA <- ggplot(as.data.frame(de.results), aes(x=logFC, y=-log10(adj.P.Val),
                                                   color=gene_type,
                                                   text= paste0("</br> Gen: " ,gene_name,
                                                                "</br> Anotacion: ",annotation))) +
  geom_point(size = 1) +
  scale_colour_manual(values = volcol) + 
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"),
        axis.line.x.bottom = element_line(color = 'black'),
        axis.line.y.left   = element_line(color = 'black'),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 19),
        axis.text=element_text(size=18),
        axis.title=element_text(size=22,face="bold")) #+ 
# ggtitle("Volcano plot")

ggsave2("Volcano.png",width = 8,height = 6)


# Interactive volcano plot:
interactive_volcano <- ggplotly(volcano_RNA,x = ~logFC, y = ~adj.P.Value,
                                tooltip = "text",
                                width = 600, height = 600)
htmlwidgets::saveWidget(as_widget(interactive_volcano), "interactive_volcano.html")







# For a specific gene analysis we design a barplot function for interesting genes.
barplot.gene <- function(gene.id,gene.name,gene.expression)
{
  expression.ll <- mean(unlist(gene.expression[gene.id,
                                               c("LL_1", "LL_2")]))
  expression.hl <- mean(unlist(gene.expression[gene.id,
                                               c("HL_1", "HL_2")]))
  means <- c(expression.ll,expression.hl)
  
  expression.sd.ll <- sd(unlist(gene.expression[gene.id,
                                                c("LL_1", "LL_2")]))
  expression.sd.hl <- sd(unlist(gene.expression[gene.id,
                                                c("HL_1", "HL_2")]))
  sds <- c(expression.sd.ll,expression.sd.hl)
  
  
  png(filename = paste(gene.name,paste0(gene.id,".png"),sep="_"),width = 300)
  par(lwd=3)
  xpos <- barplot(means,col=c("blue","firebrick2"),
                  names.arg = c("LL","HL"),las=2,cex.names = 1.5,
                  ylim=c(0,max(means+sds)*1.2),cex.axis = 1.5,lwd=3,
                  main=paste(c(gene.name, "-", gene.id),collapse=" "),
                  cex.main=2)
  arrows(x0 = xpos,y0 = means + sds, x1 = xpos, y1 = means - sds,length = 0.1,
         code = 3,angle=90,lwd=2)
  dev.off()
  return(list(means,means[2]/means[1]))
}


## Auxin

setwd('~/RNAseq/Phytohormones/Auxin')

barplot.gene(gene.id="kfl00051_0080",gene.name="TAA",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00109_0340",gene.name="YUCCA",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00384_0090",gene.name="NIT",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00124_0230",gene.name="ABP1",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00071_0010",gene.name="PINs",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00341_0080",gene.name="AUX",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00452_0080",gene.name="ROP",gene.expression=gene.expression)

## ABA

setwd('~//RNAseq/Phytohormones/ABA')

barplot.gene(gene.id="kfl00895_0020",gene.name="ZEP",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00188_0050",gene.name="NSY",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00075_0080",gene.name="ABAO",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00326_0060",gene.name="GTG",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00724_0050",gene.name="PP2C",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00097_0360",gene.name="SnRK2",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00015_0390",gene.name="AREB",gene.expression=gene.expression)

## JA

setwd('~//RNAseq/Phytohormones/JA')

barplot.gene(gene.id="kfl00081_0200",gene.name="LOX2",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00028_0030",gene.name="LOX2",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00195_0150",gene.name="AOS",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00107_0260",gene.name="ACX",gene.expression=gene.expression)

## SA

setwd('~//RNAseq/Phytohormones/SA')

barplot.gene(gene.id="kfl00737_0030",gene.name="ICS",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00086_0050",gene.name="CLU3",gene.expression=gene.expression)

## Cytokinin

setwd('~/RNAseq/Phytohormones/Cytokinin')

barplot.gene(gene.id="kfl00529_0020",gene.name="AHP",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00066_0050",gene.name="TypeA ARR",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00133_0040",gene.name="TypeB ARR",gene.expression=gene.expression)


## Carotenoids

setwd('~/RNAseq/Carotenoids/')

barplot.gene(gene.id="kfl00257_0100",gene.name="LCYB",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00536_0070",gene.name="LCYE",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00003_0600",gene.name="LCYBE",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00019_0320",gene.name="PSY",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00085_0130",gene.name="PDS",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00103_0130",gene.name="PDS",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00496_0070",gene.name="ZDS",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00028_0200",gene.name="BCHB",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00137_0040",gene.name="BCHB",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00876_0010",gene.name="BCHB",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00299_0120",gene.name="BCHB",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00515_0050",gene.name="BCHB",gene.expression=gene.expression)

barplot.gene(gene.id="kfl00214_0100",gene.name="ZEP",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00266_0090",gene.name="ZEP",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00895_0020",gene.name="ZEP",gene.expression=gene.expression)
barplot.gene(gene.id="kfl01023_0010",gene.name="ZEP",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00021_0010",gene.name="ZEP",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00073_0190",gene.name="ZEP",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00092_0060",gene.name="ZEP",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00604_0070",gene.name="VDE",gene.expression=gene.expression)


## PSI

setwd('~/RNAseq/PSI/')

barplot.gene(gene.id="kfl00031_0420",gene.name="psaP",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00036_0310",gene.name="psaN",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00039_0430",gene.name="psbR1",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00039_0440",gene.name="psbR1",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00057_0180",gene.name="psaE",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00193_0150",gene.name="psaD",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00213_0040",gene.name="psaF",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00269_0020",gene.name="psaH",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00283_0090",gene.name="psaO",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00319_0040",gene.name="psaG",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00519_0080",gene.name="psaK",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00792_0030",gene.name="psaP",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00881_0040",gene.name="psaL",gene.expression=gene.expression)
barplot.gene(gene.id="kfl01175_0030",gene.name="psaO-like",gene.expression=gene.expression)
barplot.gene(gene.id="kfl01813_0140",gene.name="psaCa",gene.expression=gene.expression)
barplot.gene(gene.id="kfl01813_0260",gene.name="ycf4a",gene.expression=gene.expression)
barplot.gene(gene.id="kfl01813_0270",gene.name="psaIa",gene.expression=gene.expression)
barplot.gene(gene.id="kfl01813_0320",gene.name="psaBa",gene.expression=gene.expression)
barplot.gene(gene.id="kfl01813_0330",gene.name="psaAa",gene.expression=gene.expression)
barplot.gene(gene.id="kfl01813_0900",gene.name="psaJ",gene.expression=gene.expression)
barplot.gene(gene.id="kfl01813_1260",gene.name="psaAb",gene.expression=gene.expression)
barplot.gene(gene.id="kfl01813_1270",gene.name="psaBb",gene.expression=gene.expression)
barplot.gene(gene.id="kfl01813_1320",gene.name="psaIb",gene.expression=gene.expression)
barplot.gene(gene.id="kfl01813_1330",gene.name="ycf4b",gene.expression=gene.expression)
barplot.gene(gene.id="kfl01813_1450",gene.name="psaCb",gene.expression=gene.expression)

## PSII

setwd('~/RNAseq/PSII/')

barplot.gene(gene.id="kfl00039_0430",gene.name="psbR1-1",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00039_0440",gene.name="psbR1-2",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00053_0120",gene.name="psb27",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00093_0070",gene.name="PsbS",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00098_0060",gene.name="psbR2",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00104_0070",gene.name="psbQ",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00112_0180",gene.name="psbTn",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00191_0030",gene.name="psbP1-2",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00238_0100",gene.name="psb28",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00239_0120",gene.name="psbP1-1",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00252_0060",gene.name="psb27",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00293_0040",gene.name="Psb29",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00307_0040",gene.name="HCF136",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00370_0100",gene.name="psbY",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00614_0030",gene.name="psbO",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00638_0030",gene.name="psbW",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00693_0040",gene.name="lpa3",gene.expression=gene.expression)
barplot.gene(gene.id="kfl01813_0300",gene.name="psbAa",gene.expression=gene.expression)
barplot.gene(gene.id="kfl01813_0520",gene.name="psbB",gene.expression=gene.expression)
barplot.gene(gene.id="kfl01813_0530",gene.name="psbT",gene.expression=gene.expression)
barplot.gene(gene.id="kfl01813_0540",gene.name="psbN",gene.expression=gene.expression)
barplot.gene(gene.id="kfl01813_0550",gene.name="psbH",gene.expression=gene.expression)
barplot.gene(gene.id="kfl01813_0580",gene.name="psbA",gene.expression=gene.expression)
barplot.gene(gene.id="kfl01813_0690",gene.name="psbI",gene.expression=gene.expression)
barplot.gene(gene.id="kfl01813_0710",gene.name="psbK",gene.expression=gene.expression)
barplot.gene(gene.id="kfl01813_0770",gene.name="psbM",gene.expression=gene.expression)
barplot.gene(gene.id="kfl01813_0870",gene.name="psbZ",gene.expression=gene.expression)
barplot.gene(gene.id="kfl01813_1130",gene.name="psbD",gene.expression=gene.expression)
barplot.gene(gene.id="kfl01813_1180",gene.name="psbE",gene.expression=gene.expression)
barplot.gene(gene.id="kfl01813_1190",gene.name="psbF",gene.expression=gene.expression)
barplot.gene(gene.id="kfl01813_1200",gene.name="psbL",gene.expression=gene.expression)
barplot.gene(gene.id="kfl01813_1210",gene.name="psbJ",gene.expression=gene.expression)
barplot.gene(gene.id="kfl01813_1290",gene.name="psbAb",gene.expression=gene.expression)

##Cytochrome b6/f

setwd('~/RNAseq/Cyt6f/')

barplot.gene(gene.id="kfl00195_0030",gene.name="petC2",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00433_0020",gene.name="Fe-S",gene.expression=gene.expression)
barplot.gene(gene.id="kfl01813_0470",gene.name="petG",gene.expression=gene.expression)
barplot.gene(gene.id="kfl01813_0560",gene.name="petB",gene.expression=gene.expression)
barplot.gene(gene.id="kfl01813_0570",gene.name="petD",gene.expression=gene.expression)
barplot.gene(gene.id="kfl01813_0800",gene.name="petN",gene.expression=gene.expression)

barplot.gene(gene.id="kfl00433_0020",gene.name="petC",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00005_0210",gene.name="petC",gene.expression=gene.expression)


## Calvin

setwd('~/RNAseq/Calvin/')

barplot.gene(gene.id="kfl00002_0370",gene.name="rbcS",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00093_0320",gene.name="rbcS",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00117_0060",gene.name="rbcS",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00126_0050",gene.name="PGK",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00038_0330",gene.name="PGK",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00185_0140",gene.name="PGK",gene.expression=gene.expression)

barplot.gene(gene.id="kfl00100_0310",gene.name="GAPDH",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00141_0280",gene.name="GAPDH",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00185_0150",gene.name="GAPDH",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00385_0100",gene.name="GAPDH",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00674_0010",gene.name="GAPDH",gene.expression=gene.expression)

barplot.gene(gene.id="kfl00094_0250",gene.name="TPI",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00074_0140",gene.name="TPI",gene.expression=gene.expression)

barplot.gene(gene.id="kfl00013_0290",gene.name="FBA",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00125_0270",gene.name="FBA",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00175_0130",gene.name="FBA",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00403_0050",gene.name="FBA",gene.expression=gene.expression)

barplot.gene(gene.id="kfl00027_0050",gene.name="FBPase",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00105_0190",gene.name="FBPase",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00481_0120",gene.name="FBPase",gene.expression=gene.expression)

barplot.gene(gene.id="kfl00100_0040",gene.name="Tkt",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00208_0140",gene.name="Tkt",gene.expression=gene.expression)

barplot.gene(gene.id="kfl00092_0210",gene.name="SEBP",gene.expression=gene.expression)

barplot.gene(gene.id="kfl00254_0130",gene.name="RPI",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00622_0030",gene.name="RPI",gene.expression=gene.expression)

barplot.gene(gene.id="kfl00390_0110",gene.name="PRK",gene.expression=gene.expression)

barplot.gene(gene.id="kfl00118_0220",gene.name="PSII repair",gene.expression=gene.expression)


barplot.gene(gene.id="kfl00002_0530",gene.name="AAA",gene.expression=gene.expression)



barplot.gene(gene.id="kfl00116_0080",gene.name="INDO",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00482_0070",gene.name="INDO",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00033_0260",gene.name="KMO",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00231_0160",gene.name="KMO",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00297_0040",gene.name="KMO",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00476_0010",gene.name="KMO",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00484_0050",gene.name="KMO",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00194_0090",gene.name="KYNU",gene.expression=gene.expression)

barplot.gene(gene.id="kfl00051_0080",gene.name="TAA",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00096_0120",gene.name="TAIR",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00061_0210",gene.name="TSB",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00713_0040",gene.name="TSA",gene.expression=gene.expression)


## Oxidative stress

setwd('~/RNAseq/Oxidative_stress/')

barplot.gene(gene.id="kfl00386_0040",gene.name="NTRC",gene.expression=gene.expression)
barplot.gene(gene.id="kfl01246_0010",gene.name="2CPA",gene.expression=gene.expression)

barplot.gene(gene.id="kfl00014_0250",gene.name="PRXQ",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00014_0230",gene.name="PRXQ",gene.expression=gene.expression)

barplot.gene(gene.id="kfl01057_0030",gene.name="CAT",gene.expression=gene.expression)

barplot.gene(gene.id="kfl00120_0050",gene.name="TRX",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00573_0030",gene.name="TRX",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00021_0420",gene.name="NTR",gene.expression=gene.expression)

barplot.gene(gene.id="kfl00434_0020",gene.name="eIF2",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00002_0530",gene.name="HSP90",gene.expression=gene.expression)

barplot.gene(gene.id="kfl00300_0060",gene.name="PAPST",gene.expression=gene.expression)

barplot.gene(gene.id="kfl00113_0150",gene.name="CPN60A",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00076_0150",gene.name="CPN60B",gene.expression=gene.expression)

barplot.gene(gene.id="kfl00104_0280",gene.name="NOB1",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00198_0260",gene.name="MDN1",gene.expression=gene.expression)

barplot.gene(gene.id="kfl00006_0030",gene.name="NMD3",gene.expression=gene.expression)

barplot.gene(gene.id="kfl00040_0230",gene.name="6PGDH",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00043_0100",gene.name="6PGDH",gene.expression=gene.expression)

barplot.gene(gene.id="kfl00184_0040",gene.name="EX",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00201_0150",gene.name="FtsH2",gene.expression=gene.expression)



# State-transition

setwd('~/RNAseq/state_transition/')

barplot.gene(gene.id="kfl00129_0080",gene.name="STN7",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00018_0210",gene.name="STN8",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00028_0490",gene.name="STN8",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00756_0030",gene.name="STN8",gene.expression=gene.expression)

barplot.gene(gene.id="kfl00753_0080",gene.name="PPH1",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00947_0050",gene.name="PBCP",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00118_0090",gene.name="PBCP",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00604_0010",gene.name="PBCP",gene.expression=gene.expression)


# CEF

setwd('~/RNAseq/CEF/')

barplot.gene(gene.id="kfl00020_0020",gene.name="PGR5",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00342_0140",gene.name="PGRL1",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00093_0070",gene.name="PsbS",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00478_0030",gene.name="LHCSR",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00115_0210",gene.name="NdhL",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00053_0440",gene.name="NdhM",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00564_0070",gene.name="NdhN",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00414_0120",gene.name="NdhO",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00186_0190",gene.name="NdhS",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00049_0330",gene.name="NdhT",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00355_0090",gene.name="NdhU",gene.expression=gene.expression)

barplot.gene(gene.id="kfl00651_0060",gene.name="PnsB1",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00020_0550",gene.name="PnsB3",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00267_0055",gene.name="PnsB4",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00335_0030",gene.name="PnsB5",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00214_0160",gene.name="LHCA5",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00319_0090",gene.name="CRR1",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00120_0020",gene.name="CRR6",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00076_0150",gene.name="CRR27",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00186_0190",gene.name="CRR31",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00057_0165",gene.name="CRR42",gene.expression=gene.expression)

barplot.gene(gene.id="kfl00009_0280",gene.name="PTOX",gene.expression=gene.expression)


# Protein folding

setwd('~/protein_folding/')

barplot.gene(gene.id="kfl00002_0530",gene.name="HSP90",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00090_0160",gene.name="CRT1a",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00387_0020",gene.name="CLPB3",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00475_0020",gene.name="CPN20",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00076_0150",gene.name="CPN60B",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00113_0150",gene.name="CPN60A",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00460_0030",gene.name="HSC70-7",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00006_0350",gene.name="GrpE",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00292_0140",gene.name="HSP60",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00120_0050",gene.name="PDIL",gene.expression=gene.expression)
barplot.gene(gene.id="kfl00573_0030",gene.name="TDX",gene.expression=gene.expression)


# PAP

setwd('~/RNAseq/PAP')

barplot.gene(gene.id="kfl00096_0240",gene.name="SAL1",gene.expression=gene.expression)
