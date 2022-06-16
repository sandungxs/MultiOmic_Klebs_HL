
###############################################################
#          UNIQUE STREPTOPHYTA GENES IN K. NITENS             #
###############################################################

# MARACAS Transcriptomic result:
gene.expression <- read.table(file = "gene_expression.tsv",header=T,as.is=T)
head(gene.expression)

# Mean expression:
gene.expression.ll <- (gene.expression$LL_1 + gene.expression$LL_2)/2
gene.expression.hl <- (gene.expression$HL_1 + gene.expression$HL_2)/2

mean.gene.expression <- data.frame(geneID=gene.expression$geneID,
                                   LL=gene.expression.ll,
                                   HL=gene.expression.hl)


# Activated genes of our study.
activated.genes <- read.table(file="activated_genes.txt",header=F,as.is=T)[[1]]
length(activated.genes)  # 677 activated genes
activated.mean.gene.expression <- subset(mean.gene.expression, geneID %in% activated.genes)
head(activated.mean.gene.expression)

# Repressed genes of our study.
repressed.genes <- read.table(file="repressed_genes.txt",header=F,as.is=T)[[1]]
length(repressed.genes) # 678 repressed genes
repressed.mean.gene.expression <- subset(mean.gene.expression, geneID %in% repressed.genes)
head(repressed.mean.gene.expression)



############################################
#              ORTHOLOGY STUDY             #
############################################


# BLAST orthology in Arabidopsis:
activated.genes.atha <- read.table(file="activated_genes_atha_table.txt", header=F,as.is=T)
head(activated.genes.atha)
activated.genes.atha.klebs <- activated.genes.atha$V2
names(activated.genes.atha.klebs) <- activated.genes.atha$V1

# BLAST orthology in Chlamydomonas:
activated.genes.cre <- read.table(file="activated_genes_cre_table.txt", header=F,as.is = T)
head(activated.genes.cre)
activated.genes.cre.klebs <- activated.genes.cre$V2
names(activated.genes.cre.klebs) <- activated.genes.cre$V1


# Adding orthologs ACTIVATED genes.
activated.mean.gene.expression[[4]] <- activated.genes.atha.klebs[activated.mean.gene.expression$geneID]
activated.mean.gene.expression[[5]] <- activated.genes.cre.klebs[activated.mean.gene.expression$geneID]
colnames(activated.mean.gene.expression)[4] <- "BestHitAtha"
colnames(activated.mean.gene.expression)[5] <- "BestHitCre"
head(activated.mean.gene.expression)
dim(activated.mean.gene.expression)

# K.nitens genes without C. reinhardtii homologous genes
na_cre <- activated.mean.gene.expression[is.na(activated.mean.gene.expression$BestHitCre),]
dim(na_cre)
head(na_cre)

# K.nitens genes without C. reinhardtii homologous genes and with A.thaliana homologous genes
na_cre_yes_atha <- na_cre[!is.na(na_cre$BestHitAtha),]
head(na_cre_yes_atha)
dim(na_cre_yes_atha)   # 70 unique Streptophyta genes

(nrow(na_cre_yes_atha)/length(activated.genes))*100  # 10.33 % of total activated genes.



# BLAST orthology in Arabidopsis:
repressed.genes.atha <- read.table(file="repressed_genes_atha_table.txt", header=F,as.is=T)
head(repressed.genes.atha)
repressed.genes.atha.klebs <- repressed.genes.atha$V2
names(repressed.genes.atha.klebs) <- repressed.genes.atha$V1

# BLAST orthology in Chlamydomonas:
repressed.genes.cre <- read.table(file="repressed_genes_cre_table.txt", header=F,as.is = T)
head(repressed.genes.cre)
repressed.genes.cre.klebs <- repressed.genes.cre$V2
names(repressed.genes.cre.klebs) <- repressed.genes.cre$V1

# Adding orthologs genes.
repressed.mean.gene.expression[[4]] <- repressed.genes.atha.klebs[repressed.mean.gene.expression$geneID]
repressed.mean.gene.expression[[5]] <- repressed.genes.cre.klebs[repressed.mean.gene.expression$geneID]
colnames(repressed.mean.gene.expression)[4] <- "BestHitAtha"
colnames(repressed.mean.gene.expression)[5] <- "BestHitCre"
head(repressed.mean.gene.expression)

# K.nitens genes without C. reinhardtii homologous genes
na_cre_rep <- repressed.mean.gene.expression[is.na(repressed.mean.gene.expression$BestHitCre),]
dim(na_cre_rep)
head(na_cre_rep)

# K.nitens genes without C. reinhardtii homologous genes and with A.thaliana homologous genes
na_cre_yes_atha_rep <- na_cre_rep[!is.na(na_cre_rep$BestHitAtha),]
head(na_cre_yes_atha_rep)
dim(na_cre_yes_atha_rep)   # 67 unique Streptophyta genes

(nrow(na_cre_yes_atha_rep)/length(repressed.genes))*100  # 9.88 % of total repressed genes.


