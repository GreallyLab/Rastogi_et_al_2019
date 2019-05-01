 # Get columns to filter for zcall table
Full_Data_Table_header <- read.table("Full_Data_Table_18_header.txt", sep="\t")
Full_Data_Table_header <- t(Full_Data_Table_header)
 # check file
head(Full_Data_Table_header, 20)
 # change row names to just be ordered numbers.
rownames(Full_Data_Table_header)<- 1:nrow(Full_Data_Table_header)

 # grab the columns wanted. 
fil_cols <- c("Name","Chr", ".GType","Position",".X",".Y")
idx_name <- grep(pattern = paste(fil_cols, collapse="|"), Full_Data_Table_header)
Full_Data_Table_header_fil_1 <- Full_Data_Table_header[idx_name,]

 # raw scores and the custom GType need to be removed. 
xraw_names <- c(".X Raw",".Y Raw", "Custom GType")
idx_raw <- grep(pattern= paste(xraw_names, collapse = "|"), Full_Data_Table_header_fil_1)
Full_Data_Table_header_fil_2 <- Full_Data_Table_header_fil_1[-idx_raw]

 # Visual check
head(Full_Data_Table_header_fil_2,12)
tail(Full_Data_Table_header_fil_2)

 # Sample Check
 # what are all of the included samples? Used theta because each sample only 
 # has one column matching it
header_Theta <- grep(".Theta", Full_Data_Table_header, value = T)
head(header_Theta)
length(header_Theta) # 22
str(header_Theta)
length(header_Theta)
header_Theta_substr <- do.call(rbind,strsplit(header_Theta, "[.]"))[,1]
header_Theta_substr
 # "A18" "A19" "A22" "A55" "A56" "A57" "A58" "A59" "A60" "A61" "A62" "B43" "B44" "B45" "B46"
 # "B47" "B48" "B49" "B50" "B51" "B52" "B53"
 # all of which are what we want

 # Paste a $ sign in front for easy awk filtering
Col_zcall <- paste("$", names(Full_Data_Table_header_fil_2), sep = "")
Col_zcall <- as.data.frame(Col_zcall)
head(Col_zcall)

write.table(Col_zcall, "Gtype_XY_columns.txt", quote=F,
            append = F, sep = "\t", row.names = F, col.names = F)

 # Now I want to filter the columns for LRR 
fil_LRR_cols <- c("Name","Chr", "Position","Log R Ratio")
idx_LRR_name <- grep(pattern = paste(fil_LRR_cols, collapse="|"), Full_Data_Table_header)
Full_Data_Table_header_fil_LRR <- Full_Data_Table_header[idx_LRR_name,]
length(Full_Data_Table_header_fil_LRR) # 25

 # Paste a $ sign in front for easy awk filtering
Col_LRR <- paste("$", names(Full_Data_Table_header_fil_LRR), sep = "")
Col_LRR <- as.data.frame(Col_LRR)
head(Col_LRR)

 # write out table
write.table(Col_LRR, "XY_LRR_columns.txt", quote=F,
            append = F, sep = "\t", row.names = F, col.names = F)
