options(stringsAsFactors = FALSE, scipen = 999)
library(data.table)

 # read in the filtered data tables 
ztab1 <- fread("../../MEGAarray2/Full_data_zcall_noMM_noNA_NCfil.txt")
ztab1 <- as.data.frame(ztab1)
ztab2 <- fread("Full_data_zcall_noMM_noNA_NCfil.txt")
ztab2 <- as.data.frame(ztab2)

 # remove the extra pos and chr columns
ztab2 <- ztab2[,-c(2:3)]

 # merge the tables
ztab <- merge(ztab1, ztab2, by="Name")
colnames(ztab)
dim(ztab)
 # [1] 1613843     339

# write out table
write.table(ztab,"Combined_zcall_noMM_noNA_NCfil.txt", append = F, 
            quote=F, row.names = F,
            col.names = T, sep = "\t")


