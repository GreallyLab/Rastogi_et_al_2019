### making a plot using the plink MDS 
library(scales)
options(stringsAsFactors = F)
m <- read.table("Plink_filtering/SplX_converted_noHH_indep_mds_2.mds", header=T)
plot(m[,4],m[,5], cex=.75,  col = rgb(0, 0, 1, 0.5))
head(m)

# highlight in red those samples are related.
m$color <- "black"
related <- c("B11", "B10", "A40", "A41", "B32", "B33", "B14", "B13",
             "A14", "B23", "B20", "A18", "A58", "A59")
m$color[grep(pattern= paste(related, collapse = "|"), x=m$IID)]<- "red"
m$color[which(x=m$IID=="A2")]<- "red"
m$color[which(x=m$IID=="A3")]<- "red"
m$color[which(x=m$IID=="A7")]<- "red"
m$color[which(x=m$IID=="A9")]<- "red"

pdf(file = "MDS_plot_col_by_related.pdf", width=6, height=4.5, family="ArialMT")
plot(m[,4],m[,5], cex=.75,  col = alpha(m$color, .5), ylab = "MDS 2", 
     xlab = "MDS 2", main = "Related individuals are randomly\nsparsed ancestry background")
legend("bottomright", legend=c("Related", "Unrelated"),
       pch=16, col=alpha(c("red", "black"),.5), cex=0.8)
dev.off()

# highlight in red those samples with high homozygosity
m$color_h <- "black"
homozyg <- c("A23", "B17", "B16", "A30", "A48")
m$color_h[grep(pattern= paste(homozyg, collapse = "|"), x=m$IID)]<- "red"
m$color_h[which(x=m$IID=="A2")]<- "red"
m$color_h[which(x=m$IID=="A1")]<- "red"
pdf(file = "MDS_plot_col_by_homozygosity.pdf", width=6, height=4.5, family="ArialMT")
plot(m[,4],m[,5], cex=.75,  col = alpha(m$color_h, .5), ylab = "MDS 2", 
     xlab = "MDS 2", main = "High Homozygosity samples due to distance\nfrom African ancestry axis")
legend("bottomright", legend=c("High Homozyg.", "Low Homozyg"),
       pch=16, col=alpha(c("red", "black"),.5), cex=0.8)
dev.off()

