#############################################################
#          SATURATION STUDY FOR KLEBS IN HL                 #
#############################################################

# Necessary packages:
library(ggplot2)
library(cowplot)

## Klebsormidium N2

gene.expressionN2 <- read.table(file = "data_klebs/gene_expression_N2.tsv",header=T,as.is=T)

gene.ids <- gene.expressionN2$geneID
gene.expressionN2 <- as.matrix(gene.expressionN2[,2:ncol(gene.expressionN2)])
rownames(gene.expressionN2) <- gene.ids
head(gene.expressionN2)

# Total gene expression percentage.
kN2<-(sum(apply(X = gene.expressionN2,MARGIN = 1,FUN = max) > 0)/nrow(gene.expressionN2))*100 # (11469 expressed genes) 66.33%

# Activated gene percentage
activate.genes.N2 <- read.table(file = "data_klebs/activated_genes_N2.txt",header=F, as.is = T)[[1]]
head(activate.genes.N2)
length(activate.genes.N2) # 489 activated genes
actkN2<-489/nrow(gene.expressionN2)*100 # 2.83 %

# Repressed gene percentage
repressed.genes.N2 <- read.table(file = "data_klebs/repressed_genes_N2.txt",header=F, as.is = T)[[1]]
head(repressed.genes.N2)
length(repressed.genes.N2) # 404 repressed genes
repkN2<-404/nrow(gene.expressionN2)*100  # 2.37 %


# Klebsormidium N4

gene.expressionN4 <- read.table(file = "data_klebs/gene_expression_N4.tsv",header=T,as.is=T)

gene.ids <- gene.expressionN4$geneID
gene.expressionN4 <- as.matrix(gene.expressionN4[,2:ncol(gene.expressionN4)])
rownames(gene.expressionN4) <- gene.ids
head(gene.expressionN4)

# Total gene expression percentage.
kN4<-(sum(apply(X = gene.expressionN4,MARGIN = 1,FUN = max) > 0)/nrow(gene.expressionN4))*100 # (11650 expressed genes) 67.38%
# (677 + 678)/nrow(gene.expressionN4)

# Activated gene percentage
activate.genes.N4 <- read.table(file = "data_klebs/activated_genes_N4.txt",header=F, as.is = T)[[1]]
length(activate.genes.N4) # 560 activated genes
actkN4 <- 560/nrow(gene.expressionN4)*100 # 3.24 %

# Repressed gene percentage
repressed.genes.N4 <- read.table(file = "data_klebs/repressed_genes_N4.txt",header=F, as.is = T)[[1]]
length(repressed.genes.N4)  # 533 repressed genes
repkN4<-533/nrow(gene.expressionN4)*100 # 3.08 %



# Klebsormidium N6

gene.expressionN6 <- read.table(file = "data_klebs/gene_expression_N6.tsv",header=T,as.is=T)

gene.ids <- gene.expressionN6$geneID
gene.expressionN6 <- as.matrix(gene.expressionN6[,2:ncol(gene.expressionN6)])
rownames(gene.expressionN6) <- gene.ids
head(gene.expressionN6)

# Total gene expression percentage.
kN6 <-(sum(apply(X = gene.expressionN6,MARGIN = 1,FUN = max) > 0)/nrow(gene.expressionN6))*100 # (11727 expressed genes) 67.83%

# Activated gene percentage
activate.genes.N6 <- read.table(file = "activated_genes_N6.txt",header=F, as.is = T)[[1]]
length(activate.genes.N6) # 612 activated genes
actkN6 <-612/nrow(gene.expressionN6)*100 # 3.54 %

# Repressed gene percentage
repressed.genes.N6 <- read.table(file = "data_klebs/repressed_genes_N6.txt",header=F, as.is = T)[[1]]
length(repressed.genes.N6) # 584 repressed genes
repkN6 <-584/nrow(gene.expressionN6)*100 # 3.38 %



# Klebsormidium N8

gene.expressionN8 <- read.table(file = "data_klebs/gene_expression_N8.tsv",header=T,as.is=T)

gene.ids <- gene.expressionN8$geneID
gene.expressionN8 <- as.matrix(gene.expressionN8[,2:ncol(gene.expressionN8)])
rownames(gene.expressionN8) <- gene.ids
head(gene.expressionN8)

# Total gene expression percentage.
kN8 <- (sum(apply(X = gene.expressionN8,MARGIN = 1,FUN = max) > 0)/nrow(gene.expressionN8))*100 # (11752 expressed genes) 67.97%

# Activated gene percentage
activate.genes.N8 <- read.table(file = "data_klebs/activated_genes_N8.txt",header=F, as.is = T)[[1]]
length(activate.genes.N8) # 634 activated genes
actkN8 <- 634/nrow(gene.expressionN8)*100 # 3.67 %

# Repressed gene percentage
repressed.genes.N8 <- read.table(file = "data_klebs/repressed_genes_N8.txt",header=F, as.is = T)[[1]]
length(repressed.genes.N8) # 621 repressed genes
repkN8 <- 621/nrow(gene.expressionN8)*100 # 3.59 %


# Klebsormidium N10

gene.expressionN10 <- read.table(file = "data_klebs/gene_expression_N10.tsv",header=T,as.is=T)

gene.ids <- gene.expressionN10$geneID
gene.expressionN10 <- as.matrix(gene.expressionN10[,2:ncol(gene.expressionN10)])
rownames(gene.expressionN10) <- gene.ids
head(gene.expressionN10)

# Total gene expression percentage.
kN10 <- (sum(apply(X = gene.expressionN10,MARGIN = 1,FUN = max) > 0)/nrow(gene.expressionN10))*100 # (11773 expressed genes) 68.09%

# Activated gene percentage
activate.genes.N10 <- read.table(file = "data_klebs/activated_genes_N10.txt",header=F, as.is = T)[[1]]
length(activate.genes.N10)  # 669 activated genes
actkN10 <- 669/nrow(gene.expressionN10)*100 # 3.87 %

# Repressed gene percentage
repressed.genes.N10 <- read.table(file = "data_klebs/repressed_genes_N10.txt",header=F, as.is = T)[[1]]
length(repressed.genes.N10)  # 649 repressed genes
repkN10 <- 649/nrow(gene.expressionN10)*100 # 3.75 %


# Klebsormidium N12

gene.expressionN12 <- read.table(file = "data_klebs/gene_expression_N12.tsv",header=T,as.is=T)

gene.ids <- gene.expressionN12$geneID
gene.expressionN12 <- as.matrix(gene.expressionN12[,2:ncol(gene.expressionN12)])
rownames(gene.expressionN12) <- gene.ids
head(gene.expressionN12)

# Total gene expression percentage.
kN12 <-(sum(apply(X = gene.expressionN12,MARGIN = 1,FUN = max) > 0)/nrow(gene.expressionN12))*100 # (11796 expressed genes) 68.22%

# Activated gene percentage
activate.genes.N12 <- read.table(file = "data_klebs/activated_genes_N12.txt",header=F, as.is = T)[[1]]
length(activate.genes.N12) # 684 activated genes
actkN12 <- 684/nrow(gene.expressionN12)*100 # 3.96 %

# Repressed gene percentage
repressed.genes.N12 <- read.table(file = "data_klebs/repressed_genes_N12.txt",header=F, as.is = T)[[1]]
length(repressed.genes.N12) # 656 repressed genes
repkN12 <- 656/nrow(gene.expressionN12)*100 # 3.79 %


# Klebsormidium N14

gene.expressionN14 <- read.table(file = "data_klebs/gene_expression_N14.tsv",header=T,as.is=T)

gene.ids <- gene.expressionN14$geneID
gene.expressionN14 <- as.matrix(gene.expressionN14[,2:ncol(gene.expressionN14)])
rownames(gene.expressionN14) <- gene.ids
head(gene.expressionN14)

# Total gene expression percentage.
kN14 <- (sum(apply(X = gene.expressionN14,MARGIN = 1,FUN = max) > 0)/nrow(gene.expressionN14))*100 # (11807 expressed genes) 68.28%

# Activated gene percentage
activate.genes.N14 <- read.table(file = "data_klebs/activated_genes_N14.txt",header=F, as.is = T)[[1]]
length(activate.genes.N14) # 672 activated genes
actkN14 <- 672/nrow(gene.expressionN14)*100 # 3.88%

# Repressed gene percentage
repressed.genes.N14 <- read.table(file = "data_klebs/repressed_genes_N14.txt",header=F, as.is = T)[[1]]
length(repressed.genes.N14) # 688 repressed genes
repkN14 <- 688/nrow(gene.expressionN14)*100 # 3.98 %



# Klebsormidium N16

gene.expressionN16 <- read.table(file = "data_klebs/gene_expression_N16.tsv",header=T,as.is=T)

gene.ids <- gene.expressionN16$geneID
gene.expressionN16 <- as.matrix(gene.expressionN16[,2:ncol(gene.expressionN16)])
rownames(gene.expressionN16) <- gene.ids
head(gene.expressionN16)

# Total gene expression percentage.
kN16 <- (sum(apply(X = gene.expressionN16,MARGIN = 1,FUN = max) > 0)/nrow(gene.expressionN16))*100 # (11814 expressed genes) 68.32%

# Activated gene percentage
activate.genes.N16 <- read.table(file = "data_klebs/activated_genes_N16.txt",header=F, as.is = T)[[1]]
length(activate.genes.N16) # 676 activated genes
actkN16 <- 676/nrow(gene.expressionN16)*100 # 3.91 %

# Repressed gene percentage
repressed.genes.N16 <- read.table(file = "data_klebs/repressed_genes_N16.txt",header=F, as.is = T)[[1]]
length(repressed.genes.N16) # 681 repressed genes
repkN16 <- 681/nrow(gene.expressionN16)*100 # 3.93%



# Klebsormidium N17

gene.expressionN17 <- read.table(file = "data_klebs/gene_expression_N17.tsv",header=T,as.is=T)

gene.ids <- gene.expressionN17$geneID
gene.expressionN17 <- as.matrix(gene.expressionN17[,2:ncol(gene.expressionN17)])
rownames(gene.expressionN17) <- gene.ids
head(gene.expressionN17)

# Total gene expression percentage.
kN17 <- (sum(apply(X = gene.expressionN17,MARGIN = 1,FUN = max) > 0)/nrow(gene.expressionN17))*100 # (11813 expressed genes) 68.32%
# (677 + 678)/nrow(gene.expressionN17)

# Activated gene percentage
activate.genes.N17 <- read.table(file = "data_klebs/activated_genes_N17.txt",header=F, as.is = T)[[1]]
length(activate.genes.N17)  # 686 activated genes
actkN17 <- 686/nrow(gene.expressionN17)*100 # 3.97 %

# Repressed gene percentage
repressed.genes.N17 <- read.table(file = "data_klebs/repressed_genes_N17.txt",header=F, as.is = T)[[1]]
length(repressed.genes.N17) # 679 repressed genes
repkN17 <- 679/nrow(gene.expressionN17)*100 # 3.93 %



# Klebsormidium Total

gene.expressionkTotal <- read.table(file = "data_klebs/gene_expression.tsv",header=T,as.is=T)

gene.ids <- gene.expressionkTotal$geneID
gene.expressionkTotal <- as.matrix(gene.expressionkTotal[,2:ncol(gene.expressionkTotal)])
rownames(gene.expressionkTotal) <- gene.ids
head(gene.expressionkTotal)

# Total gene expression percentage.
kTotal<-sum(apply(X = gene.expressionkTotal,MARGIN = 1,FUN = max) > 0)/nrow(gene.expressionkTotal)*100 # (11820 expressed genes) 68.36%

# Activated gene percentage
actkTotal <- 677/nrow(gene.expressionkTotal)*100  # 3.92 %

# Repressed gene percentage
repkTotal <- 678/nrow(gene.expressionkTotal)*100 # 3.92 %



## Plots

reads <-c(2,4,6,8,10,12,14,16,17,18)
percentages <- c(kN2,kN4,kN6,kN8,kN10,kN12,kN14,kN16,kN17,kTotal)
activated <- c(actkN2,actkN4,actkN6,actkN8,actkN10,actkN12,actkN14,actkN16,actkN17,actkTotal)
repressed <- c(repkN2,repkN4,repkN6,repkN8,repkN10,repkN12,repkN14,repkN16,repkN17,repkTotal)
df<-as.data.frame(cbind(reads,percentages,activated,repressed))
head(df)

# Necessary packages:
library(ggplot2)
library(cowplot)

# Plot for percentage of genome covered:
ggplot()+
  geom_line(data = df,mapping = aes(x =reads, y=percentages),color="purple")+
  geom_point(data = df,mapping = aes(x =reads, y=percentages),color="purple")+
ggtitle("Percentage of genome covered \n in our experiment")+
  theme(plot.title=element_text(size=15,hjust=0.5,face="bold"))+
  xlab("Million Reads")+
  ylab("Percentage")+
  theme(axis.title=element_text(size=13,hjust=0.5))

ggsave2("Genome_covered.png",width = 7,height = 8)
  
# Plot for activated and repressed genes.
ggplot()+
  geom_line(data = df,mapping = aes(x =reads, y=activated,colour="Activated"))+
  geom_point(data = df,mapping = aes(x =reads, y=activated,colour="Activated"))+
  geom_line(data = df,mapping = aes(x =reads, y=repressed,colour="Repressed"))+
  geom_point(data = df,mapping = aes(x =reads, y=repressed,colour="Repressed"))+
  ggtitle("Percentage of genes activated and \n repressed in our experiment.")+
  theme(plot.title=element_text(size=15,hjust=0.5,face="bold"))+
  scale_colour_discrete("Genes")+
  xlab("Million Reads")+
  ylab("DEGs")

ggsave2("DEGS.png",width = 7,height = 8)


  
#####################################################
#        SATURATION STUDY FOR KLEBS IN OPDA         #
#####################################################

# OPDA N4

gene.expressionOPDAN4 <- read.table(file = "art_res/gene_expressionN4.tsv",header=T,as.is=T)

gene.ids <- gene.expressionOPDAN4$geneID
gene.expressionOPDAN4 <- as.matrix(gene.expressionOPDAN4[,2:ncol(gene.expressionOPDAN4)])
rownames(gene.expressionOPDAN4) <- gene.ids
head(gene.expressionOPDAN4)

# Total gene expression percentage.
oN4<-(sum(apply(X = gene.expressionOPDAN4,MARGIN = 1,FUN = max) > 0)/nrow(gene.expressionOPDAN4))*100 # (12523 expressed genes) 72.49%

# Activated gene percentage
activate.genes.OPDAN4 <- read.table(file = "art_res/activated_genesN4.txt",header=F, as.is = T)[[1]]
head(activate.genes.OPDAN4)
length(activate.genes.OPDAN4) # 300  activated genes
actoN4<-300/nrow(gene.expressionOPDAN4)*100 # 1.74 %

# Repressed gene percentage
repressed.genes.OPDAN4 <- read.table(file = "art_res/repressed_genesN4.txt",header=F, as.is = T)[[1]]
head(repressed.genes.OPDAN4)
length(repressed.genes.OPDAN4)  # 367 repressed genes
repoN4<-367/nrow(gene.expressionOPDAN4)*100  # 2.12 %


# OPDA N8

gene.expressionOPDAN8 <- read.table(file = "art_res/gene_expressionN8.tsv",header=T,as.is=T)

gene.ids <- gene.expressionOPDAN8$geneID
gene.expressionOPDAN8 <- as.matrix(gene.expressionOPDAN8[,2:ncol(gene.expressionOPDAN8)])
rownames(gene.expressionOPDAN8) <- gene.ids
head(gene.expressionOPDAN8)

# Total gene expression percentage.
oN8<-(sum(apply(X = gene.expressionOPDAN8,MARGIN = 1,FUN = max) > 0)/nrow(gene.expressionOPDAN8))*100 # (12737 expressed genes) 73.66%

# Activated gene percentage
activate.genes.OPDAN8 <- read.table(file = "art_res/activated_genesN8.txt",header=F, as.is = T)[[1]]
head(activate.genes.OPDAN8)
length(activate.genes.OPDAN8)  # 307 activated genes
actoN8<-307/nrow(gene.expressionOPDAN8)*100 # 1.78 %

# Repressed gene percentage
repressed.genes.OPDAN8 <- read.table(file = "art_res/repressed_genesN8.txt",header=F, as.is = T)[[1]]
head(repressed.genes.OPDAN8)
length(repressed.genes.OPDAN8)  # 363 repressed genes
repoN8<-363/nrow(gene.expressionOPDAN8)*100  # 2.10 %


# OPDA N12

gene.expressionOPDAN12 <- read.table(file = "art_res/gene_expressionN12.tsv",header=T,as.is=T)

gene.ids <- gene.expressionOPDAN12$geneID
gene.expressionOPDAN12 <- as.matrix(gene.expressionOPDAN12[,2:ncol(gene.expressionOPDAN12)])
rownames(gene.expressionOPDAN12) <- gene.ids
head(gene.expressionOPDAN12)

# Total gene expression percentage.
oN12<-(sum(apply(X = gene.expressionOPDAN12,MARGIN = 1,FUN = max) > 0)/nrow(gene.expressionOPDAN12))*100 # (12861 expressed genes) 74.389%

# Activated gene percentage
activate.genes.OPDAN12 <- read.table(file = "art_res/activated_genesN12.txt",header=F, as.is = T)[[1]]
head(activate.genes.OPDAN12)
length(activate.genes.OPDAN12)  # 304 activated genes
actoN12<-304/nrow(gene.expressionOPDAN12)*100 # 1.76 %

# Repressed gene percentage
repressed.genes.OPDAN12 <- read.table(file = "art_res/repressed_genesN12.txt",header=F, as.is = T)[[1]]
head(repressed.genes.OPDAN12)
length(repressed.genes.OPDAN12) # 377 repressed genes
repoN12<-377/nrow(gene.expressionOPDAN12)*100  # 2.18 %


# OPDA N16

gene.expressionOPDAN16 <- read.table(file = "art_res/gene_expressionN16.tsv",header=T,as.is=T)

gene.ids <- gene.expressionOPDAN16$geneID
gene.expressionOPDAN16 <- as.matrix(gene.expressionOPDAN16[,2:ncol(gene.expressionOPDAN16)])
rownames(gene.expressionOPDAN16) <- gene.ids
head(gene.expressionOPDAN16)

# Total gene expression percentage.
oN16<-(sum(apply(X = gene.expressionOPDAN16,MARGIN = 1,FUN = max) > 0)/nrow(gene.expressionOPDAN16))*100 # (12861 expressed genes) 74.389%

# Activated gene percentage
activate.genes.OPDAN16 <- read.table(file = "art_res/activated_genesN16.txt",header=F, as.is = T)[[1]]
head(activate.genes.OPDAN16)
length(activate.genes.OPDAN16) # 300 activated genes
actoN16<-300/nrow(gene.expressionOPDAN16)*100 # 1.74 %

# Repressed gene percentage
repressed.genes.OPDAN16 <- read.table(file = "art_res/repressed_genesN16.txt",header=F, as.is = T)[[1]]
head(repressed.genes.OPDAN16)
length(repressed.genes.OPDAN16)  # 379 repressed genes
repoN16<-379/nrow(gene.expressionOPDAN16)*100  # 2.19 %


# OPDA N20

gene.expressionOPDAN20 <- read.table(file = "art_res/gene_expressionN20.tsv",header=T,as.is=T)

gene.ids <- gene.expressionOPDAN12$geneID
gene.expressionOPDAN20 <- as.matrix(gene.expressionOPDAN20[,2:ncol(gene.expressionOPDAN20)])
rownames(gene.expressionOPDAN20) <- gene.ids
head(gene.expressionOPDAN20)

# Total gene expression percentage.
oN20<-(sum(apply(X = gene.expressionOPDAN20,MARGIN = 1,FUN = max) > 0)/nrow(gene.expressionOPDAN20))*100 # (12871 expressed genes) 74.44%

# Activated gene percentage
activate.genes.OPDAN20 <- read.table(file = "art_res/activated_genesN20.txt",header=F, as.is = T)[[1]]
head(activate.genes.OPDAN20)
length(activate.genes.OPDAN20)  # 300 activated genes
actoN20<-300/nrow(gene.expressionOPDAN20)*100 # 1.74 %

# Repressed gene percentage
repressed.genes.OPDAN20 <- read.table(file = "art_res/repressed_genesN20.txt",header=F, as.is = T)[[1]]
head(repressed.genes.OPDAN20)
length(repressed.genes.OPDAN20)  # 379 repressed genes.
repoN20<-379/nrow(gene.expressionOPDAN20)*100  # 2.19 %



# OPDA Total

gene.expressionOPDATotal <- read.table(file = "art_res/gene_expressionTotal.tsv",header=T,as.is=T)

gene.ids <- gene.expressionOPDATotal$geneID
gene.expressionOPDATotal <- as.matrix(gene.expressionOPDATotal[,2:ncol(gene.expressionOPDATotal)])
rownames(gene.expressionOPDATotal) <- gene.ids
head(gene.expressionOPDATotal)

# Total gene expression percentage.
oNTotal<-(sum(apply(X = gene.expressionOPDATotal,MARGIN = 1,FUN = max) > 0)/nrow(gene.expressionOPDATotal))*100 # (12871 expressed genes) 74.44%

# Activated gene percentage
activate.genes.OPDATotal <- read.table(file = "art_res/activated_genesTotal.txt",header=F, as.is = T)[[1]]
head(activate.genes.OPDATotal)
length(activate.genes.OPDATotal) # 300 activated genes
actoTotal<-300/nrow(gene.expressionOPDATotal)*100 # 1.74 %

# Repressed gene percentage
repressed.genes.OPDATotal <- read.table(file = "art_res/repressed_genesTotal.txt",header=F, as.is = T)[[1]]
head(repressed.genes.OPDATotal)
length(repressed.genes.OPDATotal)  # 379 repressed genes
repoTotal<-379/nrow(gene.expressionOPDATotal)*100  # 2.19 %




## Plots

readsOPDA <-c(4,8,12,16,20,24)
percentagesOPDA <- c(oN4,oN8,oN12,oN16,oN20,oNTotal)
activatedOPDA <- c(actoN4,actoN8,actoN12,actoN16,actoN20,actoTotal)
repressedOPDA <- c(repoN4,repoN8,repoN12,repoN16,repoN20,repoTotal)
df_OPDA<-as.data.frame(cbind(readsOPDA,percentagesOPDA,activatedOPDA,repressedOPDA))
head(df_OPDA)

library(ggplot2)
library(cowplot)

# # Plot for percentage of genome covered in OPDA experiment:
ggplot()+
  geom_line(data = df_OPDA,mapping = aes(x =readsOPDA, y=percentagesOPDA),color="purple")+
  geom_point(data = df_OPDA,mapping = aes(x =readsOPDA, y=percentagesOPDA),color="purple")+
  ggtitle("Percentage of genome covered \n in OPDA experiment")+
  theme(plot.title=element_text(size=15,hjust=0.5,face="bold"))+
  xlab("Million Reads")+
  ylab("Percentage")+
  theme(axis.title=element_text(size=13,hjust=0.5))

ggsave2("Genome_covered_OPDA.png",width = 7,height = 8)

# Plot para activados y reprimidos del experimento OPDA
ggplot()+
  geom_line(data = df_OPDA,mapping = aes(x =readsOPDA, y=activatedOPDA,colour="Activated"))+
  geom_point(data = df_OPDA,mapping = aes(x =readsOPDA, y=activatedOPDA,colour="Activated"))+
  geom_line(data = df_OPDA,mapping = aes(x =readsOPDA, y=repressedOPDA,colour="Repressed"))+
  geom_point(data = df_OPDA,mapping = aes(x =readsOPDA, y=repressedOPDA,colour="Repressed"))+
  ggtitle("Percentage of genes activated and \n repressed in OPDA experiment")+
  theme(plot.title=element_text(size=15,hjust=0.5,face="bold"))+
  scale_colour_discrete("Genes")+
  xlab("Million Reads")+
  ylab("DEGs")

ggsave2("DEGS_OPDA.png",width = 7,height = 8)

