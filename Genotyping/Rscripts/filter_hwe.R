 # code to filter out those SNPs breaking HWE p < 0.000001 

 #libraries
library(data.table)
 
 #read-in the original SNP file
Data_table <- fread("Full_data_zcall_noMM_noNA_NCfil.txt")
Data_table <- as.data.frame(Data_table)

 # read-in filtered SNPs
hwe_passing_tab <- fread("Full_data_hwe_check_2.bim")
hwe_passing_tab <- as.data.frame(hwe_passing_tab)

 #extract only those variants which pass filtering
idx_pass <- which(Data_table$Name %in% hwe_passing_tab$V2)
Data_table_hwe <- Data_table[idx_pass,]

 #now do this for the first array plate
 #read-in the original SNP file
Data_table_1 <- fread("../../MEGAarray2/Full_data_zcall_noMM_noNA_NCfil.txt")
Data_table_1 <- as.data.frame(Data_table_1)

 # read-in filtered SNPs
hwe_passing_tab_1 <- fread("../../MEGAarray2/Full_data_hwe_check_2.bim")
hwe_passing_tab_1 <- as.data.frame(hwe_passing_tab_1)

 #extract only those variants which pass filtering
idx_pass_1 <- which(Data_table_1$Name %in% hwe_passing_tab_1$V2)
Data_table_1_hwe <- Data_table_1[idx_pass_1,]

 # now merge the two tables
Data_table_com <- merge(Data_table_1_hwe, Data_table_hwe, by="Name", sort=FALSE, all=FALSE)
dim(Data_table_com) # 1602747     341

 # clean up columns
idx_rmCol <- which(colnames(Data_table_com) %in% c("Chr.y","Position.y"))
Data_table_com <- Data_table_com[,-idx_rmCol]

 # write out the file
write.table(Data_table_com,"Full_data_zcall_noMM_noNA_NCfil_noHWE_mergePlates.txt", append = F, 
            quote=F, row.names = F, col.names = T, sep = "\t")

