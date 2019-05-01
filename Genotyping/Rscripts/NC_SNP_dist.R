 # Look at the number of non-calls (NCs) in the genotype data

 # set options
options(scipen = 999, stringsAsFactors = F)
library(data.table)

 # read in the NC rows
Nc_rows <- fread("Full_data_zcall_noMM_noNA_NC.txt", header=FALSE)
 #Nc_rows <- fread("../../MEGAarray2/Full_data_zcall_noMM_noNA_NC.txt", header=FALSE)

Nc_rows <- as.data.frame(Nc_rows)
nrow(Nc_rows) # 79817

 # calculate the number of NCs per 
NC_per_row <- apply(X = Nc_rows, MARGIN = 1, FUN = function(x){sum(x=="NC")})
length(NC_per_row)

pdf(file = "NC_per_SNP_dist.pdf", width = 7, height = 5, family="arialMT") 
hist(NC_per_row, main="Distribution of NCs per SNPs with an NC", xlab = "# of Samples")
 # there are a good 30k SNPs for which all samples have NC.
dev.off()

 # get the SNPs to remove (all SNPs with >1 NC)
Nc_rows$NC <- NC_per_row
rm_NC_Snps <- Nc_rows[Nc_rows$NC>1, 1]
#rm_NC_Snps <- Nc_rows[Nc_rows$NC>4, 1]

head(rm_NC_Snps)
nrow(rm_NC_Snps) # 47,727
 # length(rm_NC_Snps)

 # Read in the full table
Zcall_data_table <- fread("Full_data_zcall_noMM_noNA.txt", header=T)
 #Zcall_data_table <- fread("../../MEGAarray2/Full_data_zcall_noMM_noNA.txt", header=T)
 # 5302 snp removed 

Zcall_data_table <- as.data.frame(Zcall_data_table)

 # filter out the >1 NC SNPs
idx_NC_snp <- which(Zcall_data_table[,1] %in% rm_NC_Snps$V1)
 # idx_NC_snp <- which(Zcall_data_table[,1] %in% rm_NC_Snps)


Zcall_data_table_filtered <- Zcall_data_table[-idx_NC_snp,]
write.table(Zcall_data_table_filtered,"Full_data_zcall_noMM_noNA_NCfil.txt", append = F, 
            quote=F, row.names = F,
            col.names = T, sep = "\t")

 #write.table(Zcall_data_table_filtered,"../../MEGAarray2/Full_data_zcall_noMM_noNA_NCfil.txt", append = F, 
 #           quote=F, row.names = F,
 #           col.names = T, sep = "\t")


