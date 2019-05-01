# testing the batch effect. 

snp_batch_1 <- read.table("../../MEGAarray2/Full_data_zcall_18_test.txt")
snp_batch_2 <- read.table("Full_data_zcall_18_test.txt")
snp_batch <- merge(snp_batch_1, snp_batch_2, by="V1")
snp_batch <- snp_batch[,-c(274,275)]
dim(snp_batch)
snp_x <- snp_batch[,seq(5, 339, 3)]
snp_y <- snp_batch[,seq(6, 339, 3)]
snp_bat <- c(rep("1", 90), rep("2", 22))
snp_geno <- snp_batch[,seq(4, 339, 3)]
snp <- cbind(t(snp_x), t(snp_y),t(snp_geno),snp_bat)
head(snp)
plot(snp, col=snp_bat)


snp_batch_3 <- read.table("../../MEGAarray2/Full_data_zcall_AB_test.txt")
snp_batch_3[12,157]
header <-read.table("../../MEGAarray2/Full_data_zcall_header.txt")
header[,157]
dim(header)
head(snp_batch_3)
snp_batch_4 <- read.table("Full_data_zcall_AB_test.txt")
snp_batch <- merge(snp_batch_3, snp_batch_4, by="V1")
snp_batch <- snp_batch[,-c(274,275)]
dim(snp_batch)
snp_bat <- c(rep("1", 90), rep("2", 22))

which(t(snp_batch[4,])=="NC")
snp_batch[4,157]
length(seq(4,157, 3))
which(t(snp_batch[16,])=="NC")


i<-16
i<-4
for(i in 1:nrow(snp_batch)){
  snp_x <- snp_batch[i,seq(5, 339, 3)]
  snp_y <- snp_batch[i,seq(6, 339, 3)]
  snp_geno <- snp_batch[i,seq(4, 339, 3)]
  snp <- cbind(t(snp_x), t(snp_y),t(snp_geno),snp_bat)
  head(snp)
  plot(snp, col=snp_bat)
  text(snp, labels=snp[,3], cex= 0.7)
}

# i <-4
test <- read.table("../../MEGAarray2/remove_this.txt")
test <- test[,-c(1:4)]
head(test)
test[,103:104]
2*52


# i<-16
test <- read.table("remove_this.txt")
test <- test[,-c(1:4)]
head(test)
test[189:190]
