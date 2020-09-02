#
# Need to remake plots for Keith's paper
# Date:  14-Apr-2020
# Name:  Kevin R Brown
#

# Set working directory
setwd("~/Documents/Work/CRISPR/Keith_CTLpaper/");


# Load data
data1 <- read.csv("result_table_diffGI_data_IFN_endTP.txt",header=T,sep="\t",stringsAsFactors=F);
data1$color <- ifelse( data1$fdr < 0.01,ifelse( data1$diff < 0,"blue","yellow"),"grey");

# Start making plot
require(ggplot2);
require(ggrepel);

ggplot( data1, aes(x=WT,y=Fitm2)) +
       geom_point( aes( size=fdr,color=color ),shape=20 ) +
       scale_color_manual(values=c('blue','grey', 'yellow'))+
       scale_size( trans="reciprocal",range=c(0.1,10)) +
       geom_smooth( method=lm, linetype="solid", color="red", size=0.5 ) +
       theme(legend.title = element_text( size=12,face="bold"),
             panel.grid.major = element_line(size = 0.25, linetype = 'solid',colour = "grey"),
             panel.grid.minor = element_blank(),
             panel.background = element_rect(fill = "#FFFFFF", colour = "black",size = 1, linetype = "solid"),
             axis.text.x = element_text( color="black",size=12,angle=45,hjust=1),
             axis.text.y = element_text( color="black",size=12),
             axis.title = element_text( color="black",size=14,face="bold"),
             plot.title = element_text( size=16, face="bold", hjust=0.5 ));


THRESH <- 0.01;
ggplot( data = data1 ) +
        geom_point( aes(x=WT, y=Fitm2),shape=20,size=0.1,colour="lightgrey" ) +
        geom_point( data=subset( data1, fdr < THRESH & diff < 0 ),aes(x=WT,y=Fitm2,size=-log10( fdr )),pch=21,colour="black",fill="blue" ) + 
        geom_point( data=subset( data1, fdr < THRESH & diff > 0 ),aes(x=WT,y=Fitm2,size=-log10( fdr )),pch=21,colour="black",fill="yellow" ) +
        geom_text_repel( data=subset( data1, fdr < THRESH ),aes(x=WT,y=Fitm2,label=X),size=2) + 
        labs( title = "Fitm2 - IFN Tx - End-point" ) + 
        xlab("WT (NormZ)") +
        ylab("Fitm2 (NormZ)") +
        geom_abline( slope = 1, intercept = 0, colour="red" ) +
       theme(legend.position = "none",
             panel.grid.major = element_line(size = 0.25, linetype = 'solid',colour = "grey"),
             panel.grid.minor = element_blank(),
             panel.background = element_rect(fill = "#FFFFFF", colour = "black",size = 1, linetype = "solid"),
             axis.text.x = element_text( color="black",size=12,angle=45,hjust=1),
             axis.text.y = element_text( color="black",size=12),
             axis.title = element_text( color="black",size=14,face="bold"),
             plot.title = element_text( size=16, face="bold", hjust=0.5 ));
ggsave( "~/Dropbox/CTL_paper/Revised Manuscript/New figures/KRB/scatterplot_Fitm2_IFNtx_endPoint.pdf",plot=last_plot(),width=5,height=5,bg="transparent",dpi=300,useDingbats=F);