options(stringsAsFactors = F, scipen = 999)

pca <- read.table("genes_90percent.pca", header=T)
head(pca)

## outputting the first 3 PCs to add to covariate table
head(pca,1)
pca[,1] <- paste("Gene_exprs_", 
                       do.call(rbind,strsplit(pca[,1], split="_"))[,6],
                       sep="")
colnames(pca)[-1] <- substring(colnames(pca)[-1], 2, 7)
cov3 <- pca[1:3,]
head(cov3)
write.table(cov3, "Gene_exprs_cov_3PCs.txt", col.names = T, row.names = F,
            append = F, quote = F, sep = "\t")


