# Make the LRR_plot

 # setup
library(data.table)
library(scales)
options(scipen=999, stringsAsFactors = F)
XY_table <- fread("Full_data_XY_LRR.txt")
XY_table <- as.data.frame(XY_table)
nrow(XY_table) # 55575

 #filter out the Multimapping and NA snps
snp_pass_qc <- fread("Full_data_zcall_noMM_noNA_noNC.snps", header=T, sep="\t")
snp_pass_qc <- as.data.frame(snp_pass_qc)

idx_pass <- which(XY_table$V1 %in% snp_pass_qc$Name)
length(idx_pass) # 43983

XY_table_pass <- XY_table[idx_pass,]

X_LRR <- XY_table_pass[XY_table_pass$V2=="X",]
dim(X_LRR)
Y_LRR <- XY_table_pass[XY_table_pass$V2=="Y",]
 # Y_LRR <- replace(Y_LRR,is.na(Y_LRR), 0) # using 0 will skew
Y_LRR<- na.omit(Y_LRR)
dim(Y_LRR) # 88

X_LRR_mean <- apply(X_LRR[,-c(1:3)], 2, function(x) {mean(x)} )
Y_LRR_mean <- apply(Y_LRR[,-c(1:3)], 2, function(x) {mean(x)} )

plot(x=X_LRR_mean, y=Y_LRR_mean)
 # appears to be a good mix of male and female.

# get sample information
samp_info <- read.csv("plate2_sex_table.csv", header=F)
head(samp_info)

samp_XY_LRR <- cbind(samp_info, X_LRR_mean, Y_LRR_mean)
head(samp_XY_LRR)
samp_XY_LRR$colors <- "blue"
samp_XY_LRR$colors[samp_XY_LRR$V2 == 1] <- "red"

pdf(file = "LRR_plot_plate2.pdf", width = 7, height = 5, family="ArialMT") 
plot(x=samp_XY_LRR$X_LRR_mean, y=samp_XY_LRR$Y_LRR_mean, 
     col = alpha(samp_XY_LRR$colors, .5), pch=16, cex=.75, xlab="Mean X intensity",
     ylab = "Mean Y intensity", main= "Samples match designated sex: Plate #2")
dev.off()
