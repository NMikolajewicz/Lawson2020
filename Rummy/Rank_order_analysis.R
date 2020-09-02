setwd("C:/Users/rumi/Desktop/")


## reading in coreCTLs
#Cores_daisy <- read.table("Screen_recent/Screens/drugZ_V1/Output/core_drugZV1.txt", header = T, stringsAsFactors = F)

## selecting mid timepoint data for all screen, except OVA_QR
col_index=which(!grepl("OVA_QR",colnames(gene_data_table.combined)) & 
                  !grepl("late",colnames(gene_data_table.combined)) &
                  !grepl("early",colnames(gene_data_table.combined)))
colnames(gene_data_table.combined)[col_index]


gene_data_table.combined_mid=gene_data_table.combined[,col_index]
rownames(gene_data_table.combined_mid)=gene_data_table.combined_mid$GENE
nrow(gene_data_table.combined_mid) 
nrow(gene_data_table.combined)  

## extracting coreCTL genes
library(dplyr)  
cores_daisy_table.combined_mid <- inner_join(Cores_daisy, gene_data_table.combined_mid, by = "GENE")
nrow(cores_daisy_table.combined_mid)



##=================================
## calculate rank-product statistics for cores_daisy_table.combined_mid
##=================================
colnames(cores_daisy_table.combined_mid)
cores_daisy_table.combined_mid.normZ=cores_daisy_table.combined_mid[,grep("normZ",colnames(cores_daisy_table.combined_mid))]
cores_daisy_table.combined_mid.rank_synth=cores_daisy_table.combined_mid[,grep("rank_synth",colnames(cores_daisy_table.combined_mid))]
cores_daisy_table.combined_mid.rank_supp=cores_daisy_table.combined_mid[,grep("rank_supp",colnames(cores_daisy_table.combined_mid))]
cores_daisy_table.combined_mid.fdr_synth=cores_daisy_table.combined_mid[,grep("fdr_synth",colnames(cores_daisy_table.combined_mid))]
cores_daisy_table.combined_mid.fdr_supp=cores_daisy_table.combined_mid[,grep("fdr_supp",colnames(cores_daisy_table.combined_mid))]



calculate_rank_product_stat=function(rank_matrix){
  rank_matrix.scaled=apply(rank_matrix,2,rank)
  rank_matrix.dim=dim(rank_matrix)
  rank_product=apply(rank_matrix.scaled,1,function(x) prod(x)^(1/rank_matrix.dim[2]))
  
  rank.random=lapply(1:(rank_matrix.dim[2]*1000),function(x) {
    set.seed(x)
    sample(1:rank_matrix.dim[1])
  })
  
  rank.random.matrix=do.call(cbind,rank.random)
  
  d=1:(rank_matrix.dim[2]*1000)
  shuffle_index=split(d, ceiling(seq_along(d)/ncol(rank_matrix)))
  rank.random.matrix=lapply(shuffle_index,function(x) rank.random.matrix[,x])
  rank.random.matrix.rp=lapply(rank.random.matrix,function(x) apply(x,1,function(y) prod(y)^(1/rank_matrix.dim[2])))
  rank.random.matrix.rp=do.call(cbind,rank.random.matrix.rp)
  rank_product.matrix=matrix(rep(rank_product,1000),ncol=1000,byrow = F)
  
  ERP=apply(rank.random.matrix.rp-rank_product.matrix<0,1,sum)/1000
  PFP=ERP/rank(rank_product)
  list(rank_product=rank_product,EPR=ERP,PFP=PFP)
}



## calculate rank statistics separately for synth and supp
rank_synth_stat=calculate_rank_product_stat(cores_daisy_table.combined_mid.rank_synth)
rank_supp_stat=calculate_rank_product_stat(cores_daisy_table.combined_mid.rank_supp)



cores_daisy_table.combined_mid.summary_table=data.frame(GENE=cores_daisy_table.combined_mid$GENE,
                                                       normZ.mean=round(apply(cores_daisy_table.combined_mid.normZ,1,function(x) mean(x,na.rm = T)),2),
                                                       normZ.min=apply(cores_daisy_table.combined_mid.normZ,1,function(x) min(x,na.rm = T)),
                                                       normZ.max=apply(cores_daisy_table.combined_mid.normZ,1,function(x) max(x,na.rm = T)),
                                                       normZ.sd=round(apply(cores_daisy_table.combined_mid.normZ,1,function(x) sd(x,na.rm = T)),2),
                                                       rank_synth.rank_product=rank(rank_synth_stat[[1]]),
                                                       rank_synth.pvalue=rank_synth_stat[[2]],
                                                       rank_supp.rank_product=rank(rank_supp_stat[[1]]),
                                                       rank_supp.pvalue=rank_supp_stat[[2]] ,
                                                       synth_sign_exps=apply(cores_daisy_table.combined_mid.fdr_synth<=0.05,1,function(x) sum(x,na.rm=T)),
                                                       supp_sign_exps=apply(cores_daisy_table.combined_mid.fdr_supp<=0.05,1,function(x) sum(x,na.rm=T)),
                                                       cores_daisy_table.combined_mid.normZ,check.names = F,stringsAsFactors = F)



###pre-processing for Plot: log(Rank_pVal) vs Rank_Product

synth <- cores_daisy_table.combined_mid.summary_table[, c("GENE","normZ.mean","rank_synth.rank_product","rank_synth.pvalue")]
colnames(synth)[3:4] <- c("Rank_Product", "Rank_pValue")
synth$LogPVal <- -log10(synth$Rank_pValue +0.001)

supp <- cores_daisy_table.combined_mid.summary_table[, c("GENE","normZ.mean","rank_supp.rank_product","rank_supp.pvalue")]
colnames(supp)[3:4] <- c("Rank_Product", "Rank_pValue")
supp$Rank_Product <- rank(-supp$Rank_Product)
supp$LogPVal <- -log10(supp$Rank_pValue +0.001)

rank_combined <- rbind(synth, supp)


##shortening off y-axis
rank_comb_yaxisCut <- subset(rank_combined, rank_combined$LogPVal > 0.5)

##significant genes
pval_sig <- subset(rank_comb_yaxisCut, rank_comb_yaxisCut$Rank_pValue < 0.05)


library(ggplot2)
library(ggrepel)

(b <- ggplot(rank_comb_yaxisCut, aes(Rank_Product, LogPVal)) + geom_point(col = "grey")+
    geom_hline(yintercept=-log10(0.05 +0.001), linetype= 4, color = "red", size=0.5)+
    geom_point(data = pval_sig, aes(Rank_Product, LogPVal, size = LogPVal), pch = 21, col = "black",
               fill = ifelse(pval_sig$Rank_Product > 100, "yellow", "blue"))+
    geom_text_repel(data = pval_sig, aes(x = Rank_Product, y = LogPVal), label = pval_sig$GENE)+
    ylab("-log10(Rank_pVal)")+
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.line = element_line(colour = "Black"),
          panel.background = element_rect(fill = "transparent", colour = NA),
          plot.background = element_rect(fill = "transparent", colour = NA)) +
    annotate(geom="text", x=100, y=-log10(0.04 +0.001), label="p.Value = 0.05",
             color="black"))




ggsave("Rank_overlap_2_Core_rank_aggregate.pdf", plot = last_plot(),
       #width = 15, height = 7,
       width = 7, height = 7,
       bg = "transparent", dpi = 300, useDingbats=FALSE)




