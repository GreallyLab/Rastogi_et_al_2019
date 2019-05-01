dataFile <- commandArgs(6)

# load in libraries
library(ggplot2)
options(stringsAsFactors = F)

# read in match file
cat(paste("Reading in: ", dataFile,".\n",sep=""))
mat <- read.table(dataFile, header=T)

# compute percent consistent for het and hom
mat$perc_het_consistent <- mat$n_het_consistent/mat$n_het_covered
mat$perc_hom_consistent <- mat$n_hom_consistent/mat$n_hom_covered

# Get the sample ID
name <- substr(dataFile, 1, 3)

# match the bam ID to the VCF ID
idx_name <- grep(name, mat$SampleID)

if (length(idx_name)==0){
  cat(paste("There is no sample in the VCF matching this ID: ", name,".\n",sep=""))
}

# Create color factor and size
mat$color <- paste("not", name, sep = " ")
mat$color[idx_name] <- name
mat$color <- factor(mat$color, levels=c(paste("not", name, sep = " "),name))

# plot the matching
pdf(file = paste(name,"_matchPlot", ".pdf",sep = ""), width=6, height=4.5, family="ArialMT")
ggplot(mat, aes(x=perc_het_consistent, y=perc_hom_consistent, colour = color))+
         geom_point(aes(alpha = color)) +
         scale_alpha_discrete(range = c(.2, 1), guide=FALSE) +
         scale_colour_manual(values=c("black", "red"),
                             name="") +
         ylab("% hom consistent") +
         xlab("% het consistent") +
         ggtitle(paste(name, "RNA/Geno match", sep=" "))
dev.off()
cat(paste(paste(name,"_matchPlot", ".pdf",sep = "")," was successfully generated.\n",sep=""))