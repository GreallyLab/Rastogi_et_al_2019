 # Looking at the concordance of ZCall and Genome Studio using the calibrate thresholds
 # load options
options(scipen = 999, stringsAsFactors = FALSE)
concord_list <- data.frame() 
for (i in 3:15) {
  temp_list <- read.table(paste("Combined_zcall_noMM_noNA_NCfil.concordance.stats.",
                                i,".txt", sep = ""),  sep=":", skip =3)
  concord_list<- rbind(concord_list,temp_list[1,])
}
concord_list <- cbind(seq(3,15),concord_list)
concord_list
concord_list[which.max(concord_list$V2),1] # 6

