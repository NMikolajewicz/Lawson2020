# Date: Sep 28, 2018
# Author: Shirley Hui
# Takes pathway themes manually identified via a Cytoscape network created using the the core killing gprofiler results (see CoreKillingGProfiler.R).  Pathways were grouped together to form themes if they contained 30% or more similar genes.  Plot the pathway themes into bar plot.
# Input: core killing gprofiler results, core killing enrichment themes
# Output: Bar plot of core killing enriched themes  
gprofilerResults<- read.delim("/Users/Keith/Desktop/Revision Docs/Pathway analysis/Output/gprofiler_results/coreCTLgenes_gprofiler.txt")
emThemes <- read.delim("/Users/Keith/Desktop/Revision Docs/Pathway analysis/Output/enr_file/coreCTLgenes_sim0.7_enrTheme.txt",header=FALSE)
themes <- as.character(unique(emThemes[,1]))
results <- c()
for (ixx in 1:length(themes)) {
   ix <- which(emThemes[,1]==themes[ixx])
   ixs <- c()
   for (i in 1:length(ix)) {
      goid <- as.character(emThemes[ix[i],2])
      ixi <- which(gprofilerResults$term.id==goid)
      ixs <- c(ixs,ixi)
   }
   mean_overlap <- mean(gprofilerResults[ixs,]$overlap.size/gprofilerResults[ixs,]$term.size)*100
   mean_overlap.size <- mean(gprofilerResults[ixs,]$overlap.size)
   mean_term.size <- mean(gprofilerResults[ixs,]$term.size)
   min_pvalue <- -log(min(gprofilerResults[ixs,]$p.value))
   results <- rbind(results,c(mean_overlap,min_pvalue,mean_overlap.size,mean_term.size))
}
rownames(results) <- themes
colnames(results) <- c("overlap","p.value","overlap.size","term.size")

library(ggplot2)
library(RColorBrewer)

#cbPalette <- c("#FED976", "#FD8D3C", "#FC4E2A", "#E31A1C", "#aa0022")
cbPalette <- c("#ededed", "#cccccc", "#969696", "#636363", "#252525") #http://colorbrewer2.org/#type=sequential&scheme=Greys&n=5
cols = cbPalette #<- brewer.pal(6, "YlOrRd")
df = data.frame(results)
df$ratio <- paste(round(df$overlap.size,1), round(df$term.size, 1), sep = "/") #this line adds the overlap/term size ratio, rounds up the term size to xx position after comma
g = ggplot(df, aes(x = reorder (rownames(df),p.value), y = overlap)) +
   ylab("Mean Percentage Overlap") +
   theme(plot.title = element_text(hjust = -0.9)) +
   geom_col(aes(fill = p.value)) +
   geom_text(aes(label = df$ratio, hjust=0))+
   scale_fill_gradientn("-log p", colours = cols, limits=c(min(df$p.value), max(df$p.value))) +
   scale_y_continuous(position = "right") +
   theme(panel.background = element_rect(fill = "white"), axis.line.x = element_line(color="black"), axis.line.y = element_line(color="white"), axis.title.y = element_blank(),axis.ticks.y = element_blank()) +
   coord_flip()

# Adjust aspect_ratio, height and width to output figure to pdf in the desired dimensions
aspect_ratio <- 1.75
ggsave(file="/Users/Keith/Desktop/Revision Docs/Pathway analysis/Output/coreCTLgenes_sim0.7.pdf",g, height = 4 , width = 5 * aspect_ratio)





