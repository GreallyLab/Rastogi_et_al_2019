 #fixing the tfam with appropriate sample data

 # read in the zcall generated tfam file
tfam <- read.table("Merged_plates_final.tfam")
head(tfam)
nrow(tfam) #112
 # read in sample information from deepa
fam_deepa <- read.csv("Sample_table_forPLINK_all.csv")
head(fam_deepa)

 # get correct order of samples for fam file
idx_new <- match(tfam$V1, fam_deepa$Sample.ID)
fam_deepa_2 <- fam_deepa[idx_new,]
head(fam_deepa_2)
dim(fam_deepa_2)
fam_deepa_2
 # completely seperated by sex except for A7 (which is a mismatched sample)

 # switch the sample.ID column with individual ID and write out
write.table(fam_deepa_2[,c(2,1,4:7)], "Merged_plates_final.tfam", quote=F,
            append = F, sep = "\t", row.names = F, col.names = F)
