#
# Need to remake plots for Keith's paper
# Date:  14-Apr-2020
# Name:  Kevin R Brown
#

# Set working directory
setwd("~/Documents/Work/CRISPR/Keith_CTLpaper/");

# Z-scoring function
zscore <- function( x ) {
   m <- mean(x,na.rm=T);
   s <- sd(x,na.rm=T);
   (x - m)/s;
}

# Load data
data1 <- read.csv("result_table_Atg12_diffGI_data_CTL_endTP.txt",header=T,sep="\t",stringsAsFactors=F)
data2 <- read.csv("result_table_Atg12_diffGI_data_TNF_endTP.txt",header=T,sep="\t",stringsAsFactors=F);

m <- merge( data1, data2,by.x="X",by.y="X" );
m$myDiff <- m$z.y + m$z.x

# Try fisher's method on the p-values?!
tmp <- fisher.method( m[,c("pvalue.x","pvalue.y")],p.corr="none")
tmp$p.adj <- p.adjust( tmp$p.value,method="fdr" );
tmp$p.adj[ tmp$p.adj == 0 ] <- min( tmp$p.adj[ tmp$p.adj > 0 ] );
m$fdr <- tmp$p.adj;

# Start making plot
require(ggplot2);
require(ggrepel);


THRESH <- 0.01;
ggplot( data = m ) +
        geom_point( aes(x=z.x, y= z.y),shape=20,size=0.1,colour="lightgrey" ) +
        scale_x_continuous( breaks=c(-10,-5,0,5,10),limits=c(-11,12) ) + 
        scale_size(range=c(0.1,5),trans="identity") +
        geom_point( data=subset( m, fdr < THRESH & myDiff < 0 ),aes(x=z.x,y= z.y,size=-log10( fdr )),pch=21,colour="black",fill="blue" ) +
        geom_point( data=subset( m, fdr < THRESH & myDiff > 0 ),aes(x=z.x,y= z.y,size=-log10( fdr )),pch=21,colour="black",fill="yellow" ) +
        geom_text_repel( data=subset( m, fdr < 0.01 & (myDiff > 6.5 | myDiff < -7)),aes(x=z.x,y= z.y,label=X),size=3,segment.size=0.3) + 
        labs( title = "Atg12 - End-point" ) + 
        xlab("CTL_Atg12 (z, diff GI)") +
        ylab("TNF_Atg12 (z, diff GI)") +
#        geom_abline( slope = -1, intercept = 0, colour="red" ) +
       theme(legend.position = "none",
#             panel.grid.major = element_line(size = 0.25, linetype = 'solid',colour = "grey"),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             panel.background = element_rect(fill = "#FFFFFF", colour = "black",size = 1, linetype = "solid"),
             axis.text = element_text( color="black",size=12),
             axis.title = element_text( color="black",size=14,face="bold"),
             plot.title = element_text( size=16, face="bold", hjust=0.5 ));
ggsave( "~/Dropbox/CTL_paper/Revised Manuscript/New figures/KRB/scatterplot_Atg12_CTLvsTNFtx_endPoint.pdf",plot=last_plot(),width=5,height=5,bg="transparent",dpi=300,useDingbats=F);