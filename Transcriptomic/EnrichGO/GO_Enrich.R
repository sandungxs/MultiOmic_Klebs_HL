################################################
##      RNASEQ GO TERM ENRICH ANALYSIS         #
################################################

##  Necessary packages:
library(clusterProfiler,logical.return = T)
library(org.At.tair.db,logical.return = T)


## Transcriptomic data from MARACAs:
gene.expression <- read.table(file = "gene_expression.tsv",header=T,as.is=T)
head(gene.expression)
gene.ids <- gene.expression$geneID
gene.expression <- as.data.frame(gene.expression[,2:ncol(gene.expression)])
rownames(gene.expression) <- gene.ids
head(gene.expression)

gene.expression.ll <- (gene.expression$LL_1 + gene.expression$LL_2)/2
gene.expression.hl <- (gene.expression$HL_1 + gene.expression$HL_2)/2

mean.gene.expression <- data.frame(geneID=gene.ids,
                                   LL=gene.expression.ll,
                                   HL=gene.expression.hl)


############################################
#      GO Enrich for Activated Genes       #
############################################

# Activated genes:
activated.genes <- read.table(file="Activated_Enrich/activated_genes.txt",header=F,as.is=T)[[1]]
length(activated.genes)

activated.mean.gene.expression <- subset(mean.gene.expression, geneID %in% activated.genes)
head(activated.mean.gene.expression)

# BLAST orthology in Arabidopsis:
activated.genes.atha <- read.table(file="Activated_Enrich/activated_genes_atha_table.txt", header=F,as.is=T)
head(activated.genes.atha)
activated.genes.atha.klebs <- activated.genes.atha$V2
names(activated.genes.atha.klebs) <- activated.genes.atha$V1

# BLAST orthology in Chlamydomonas:
activated.genes.cre <- read.table(file="Activated_Enrich/activated_genes_cre_table.txt", header=F,as.is = T)
head(activated.genes.cre)
activated.genes.cre.klebs <- activated.genes.cre$V2
names(activated.genes.cre.klebs) <- activated.genes.cre$V1

# Adding orthologs genes.
activated.mean.gene.expression[[4]] <- activated.genes.atha.klebs[activated.mean.gene.expression$geneID]
activated.mean.gene.expression[[5]] <- activated.genes.cre.klebs[activated.mean.gene.expression$geneID]
colnames(activated.mean.gene.expression)[4] <- "BestHitAtha"
colnames(activated.mean.gene.expression)[5] <- "BestHitCre"
head(activated.mean.gene.expression)

write.table(x = activated.mean.gene.expression,file = "Activated_Enrich/activated_genes_table.tsv",quote = F,sep = "\t",row.names = F)



# GO Enrich (biological process) based on Arabidopsis universe.
atha.universe <- unique(read.table(file = "klebso2atha.tsv",header = F,sep = "\t")[[2]]) # obteined with BLAST
atha.universe <- AnnotationDbi::select(org.At.tair.db,columns = c("TAIR"),keys=keys(org.At.tair.db,keytype = "TAIR"))[[1]]

activated.atha.enrich.go <- enrichGO(gene          = activated.genes.atha$V2,
                                     #                                     universe      = atha.universe, 
                                     OrgDb         = org.At.tair.db,
                                     ont           = "BP",
                                     pAdjustMethod = "BH",
                                     pvalueCutoff  = 0.05,
                                     readable      = FALSE,
                                     keyType = "TAIR")

activated.atha.enrich.go <- as.data.frame(activated.atha.enrich.go)
head(activated.atha.enrich.go)

png(filename = "Activated_Enrich/activated_BP.png",width = 500,height = 500,res = 95)
emapplot(activated.atha.enrich.go,showCategory = 20)
dev.off()

write.table(x = activated.atha.enrich.go, file = "Activated_Enrich/activated_atha_enrich_go_universe.tsv",quote = F,sep = "\t", row.names = F)

# GO Enrich (celular component) based on Arabidopsis universe.
activated.atha.enrich.go <- enrichGO(gene          = activated.genes.atha$V2,
                                     #                                     universe      = atha.universe, 
                                     OrgDb         = org.At.tair.db,
                                     ont           = "CC",
                                     pAdjustMethod = "BH",
                                     pvalueCutoff  = 0.05,
                                     readable      = FALSE,
                                     keyType = "TAIR")


png(filename = "Activated_Enrich/activated_cc.png",width = 500,height = 500,res = 95)
emapplot(activated.atha.enrich.go)
dev.off()


activated.atha.enrich.go <- as.data.frame(activated.atha.enrich.go)
head(activated.atha.enrich.go)
write.table(x = activated.atha.enrich.go, file = "Activated_Enrich/activated_atha_enrich_go_universe_cellular_component.tsv",quote = F,sep = "\t", row.names = F)


############################################
#      GO Enrich for Repressed Genes       #
############################################

# Repressed genes:
repressed.genes <- read.table(file="Repressed_Enrich/repressed_genes.txt",header=F,as.is=T)[[1]]
length(repressed.genes)

repressed.mean.gene.expression <- subset(mean.gene.expression, geneID %in% repressed.genes)
head(repressed.mean.gene.expression)


# BLAST orthology in Arabidopsis:
repressed.genes.atha <- read.table(file="Repressed_Enrich/repressed_genes_atha_table.txt", header=F,as.is=T)
head(repressed.genes.atha)
repressed.genes.atha.klebs <- repressed.genes.atha$V2
names(repressed.genes.atha.klebs) <- repressed.genes.atha$V1

# BLAST orthology in Chlamydomonas:
repressed.genes.cre <- read.table(file="Repressed_Enrich/repressed_genes_cre_table.txt", header=F,as.is = T)
head(repressed.genes.cre)
repressed.genes.cre.klebs <- repressed.genes.cre$V2
names(repressed.genes.cre.klebs) <- repressed.genes.cre$V1

# Adding orthologs genes.
repressed.mean.gene.expression[[4]] <- repressed.genes.atha.klebs[repressed.mean.gene.expression$geneID]
repressed.mean.gene.expression[[5]] <- repressed.genes.cre.klebs[repressed.mean.gene.expression$geneID]
colnames(repressed.mean.gene.expression)[4] <- "BestHitAtha"
colnames(repressed.mean.gene.expression)[5] <- "BestHitCre"
head(repressed.mean.gene.expression)

write.table(x = repressed.mean.gene.expression,file = "Repressed_Enrich/repressed_genes_table.tsv",quote = F,sep = "\t",row.names = F)


# GO Enrich (biological process) based on Arabidopsis universe.
repressed.atha.enrich.go <- enrichGO(gene          = repressed.genes.atha$V2,
                                     #                                     universe      = atha.universe, 
                                     OrgDb         = org.At.tair.db,
                                     ont           = "BP",
                                     pAdjustMethod = "BH",
                                     pvalueCutoff  = 0.05,
                                     readable      = FALSE,
                                     keyType = "TAIR")

png(filename = "Repressed_Enrich/repressed_bp.png",width = 500,height = 500,res = 95)
emapplot(repressed.atha.enrich.go,showCategory = 20)
dev.off()

repressed.atha.enrich.go <- as.data.frame(repressed.atha.enrich.go)
head(repressed.atha.enrich.go)

write.table(x = repressed.atha.enrich.go, file = "Repressed_Enrich/repressed_atha_enrich_go_universe.tsv",quote = F,sep = "\t", row.names = F)


# GO Enrich (celular component) based on Arabidopsis universe.
repressed.atha.enrich.go <- enrichGO(gene          = repressed.genes.atha$V2,
                                     #                                     universe      = atha.universe, 
                                     OrgDb         = org.At.tair.db,
                                     ont           = "CC",
                                     pAdjustMethod = "BH",
                                     pvalueCutoff  = 0.05,
                                     readable      = FALSE,
                                     keyType = "TAIR")

png(filename = "Repressed_Enrich/repressed_cc.png",width = 500,height = 500,res = 95)
emapplot(repressed.atha.enrich.go,showCategory = 20)
dev.off()

repressed.atha.enrich.go <- as.data.frame(repressed.atha.enrich.go)
head(repressed.atha.enrich.go)

write.table(x = repressed.atha.enrich.go, file = "Repressed_Enrich/repressed_atha_enrich_go_universe_cellular_component.tsv",quote = F,sep = "\t", row.names = F)

