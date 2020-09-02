#
# Process RNASeq data for Keith's mouse project.  Alignments done
# on SciNet using STAR
#

# Set working directory
setwd("/Volumes/moffat4/kbrown/RNASeq/CCBR-RNASeq-20190716_KeithCelsius/Moffat/")

# Load data
data <- read.csv("Keith_18July19_merged_STAR_readCounts.txt",header=T,sep="\t",stringsAsFactors=F);
data <- data[ !apply( data[,-1:-2] == 0,1,all), ];

# Load annotations
annot <- read.csv("keith_samples.txt",header=T,sep="\t",stringsAsFactors=F);

# Load required libraries
require(limma);
require(edgeR);
require(rafalib);

# Look at bulk data with PCA
dge <- DGEList( counts=data[,-1:-2],genes=data[,1:2] );
grps <- cbind(annot$Cell.Line,gsub(" [ABC]$","",annot$Treatment1),gsub(" [ABC]$","",annot$Treatment2));
grps2 <- apply(grps,1,paste,collapse="_");

keep.exprs <- filterByExpr( dge,group=grps2 );
dge <- dge[keep.exprs,,keep.lib.sizes=FALSE];
dge <- calcNormFactors( dge,method="TMM" );

lcpm <- cpm( dge,log=T );
myData <- data.frame( dge$genes, lcpm );

write.table(myData,"Keith_18July19_merged_STAR_log2CPM.txt",sep="\t",quote=F,row.names=F);
pca <- prcomp( t( lcpm ) );

pdf("pcaPlot_keith_RNASeq_18July19.pdf");
cl <- as.factor( annot$Cell.Line );
bigpar();
plot( pca$x[,1],pca$x[,2],pch=18,col=rainbow(4)[as.numeric(cl)],xlab="PC1",ylab="PC2",cex=1.5,main="PCA - All Samples" );
text(c(-165,-165,75),c(115,-100,-20),labels=c("B16F10-Ova-Cas9","MC38-Ova-Cas9","Renca/HA"),cex=1.5)
dev.off();

# Comparisons:
#  1. Intergenic control [2, 3, 4] vs. Atg12-1 control [8, 9, 10]
dge <- DGEList( counts=data[,c("Moffat_01","Moffat_02","Moffat_03","Moffat_07","Moffat_08","Moffat_09")], genes=data[,1:2] );
g <- rep(c("MC38.OVA_Ctrl","MC38.OVA_Atg12.1_control"),each=3);
f <- factor( g,levels=c("MC38.OVA_Ctrl","MC38.OVA_Atg12.1_control") );

keep.exprs <- filterByExpr( dge, group=f );
dge <- dge[keep.exprs,,keep.lib.sizes=FALSE];
dge <- calcNormFactors( dge,method="TMM" );

design <- model.matrix( ~ f );

y <- voom( dge, design );
fit <- lmFit( y, design );
fit2 <- eBayes( fit );

ttable <- topTable( fit2, coef="fMC38.OVA_Atg12.1_control",number=Inf );
write.table( ttable,"topTable_Comp1_MC38.OVA_Atg12.1_Ctrl.txt",sep="\t",row.names=F );

dt <- decideTests( fit2 );
pdf("scatterplot_mc38Ova_Atg12.1_Ctrl_18July19.pdf");
bigpar();
plotMD( fit2,column=2,status=dt[,2],main="MC38OVA - Atg12.1 Control");
with(subset(ttable,adj.P.Val <= 0.05),text(AveExpr,logFC,labels=EntrezGene.ID,pos=4,cex=0.7));
abline(h=0,col="darkgrey");
dev.off();





#  Intergenic control [2, 3, 4] vs. Intergenic + IFNγ [5, 6, 7]
dge <- DGEList( counts=data[,c("Moffat_01","Moffat_02","Moffat_03","Moffat_04","Moffat_05","Moffat_06")], genes=data[,1:2] );
g <- rep(c("MC38.OVA_Ctrl","MC38.OVA_IFNg"),each=3);
f <- factor( g );

keep.exprs <- filterByExpr( dge, group=f );
dge <- dge[keep.exprs,,keep.lib.sizes=FALSE];
dge <- calcNormFactors( dge,method="TMM" );

design <- model.matrix( ~ f );

y <- voom( dge, design );
fit <- lmFit( y, design );
fit2 <- eBayes( fit );

ttable <- topTable( fit2, coef="fMC38.OVA_IFNg",number=Inf );
write.table( ttable,"topTable_Comp2_MC38.OVA_IFNg_vs_Ctrl.txt",sep="\t",row.names=F );

dt <- decideTests( fit2 );
pdf("scatterplot_mc38Ova_intergenic_IFNg_Tx_18July19.pdf");
bigpar();
plotMD( fit2,column=2,status=dt[,2],main="MC38OVA - Intergenic sgRNA - +/- IFNg");
with(subset(ttable,adj.P.Val <= 0.05),text(AveExpr,logFC,labels=EntrezGene.ID,pos=4,cex=0.7));
abline(h=0,col="darkgrey");
dev.off();



#  3. Atg12-1 control [8, 9, 10] vs. Atg12-1 + IFNγ [11, 12, 13]
#      The PCA shows that Moffat_12, which is supposed to be MC38-Ova, clusters
#      with the Renca-HA samples, while the Moffat_15 samples (identical in PCA) 
#      cluster with MC38-Ova.  Therefore, we're assuming those samples were flipped,
#      and doing the analysis accordingly
dge <- DGEList( counts=data[,c("Moffat_07","Moffat_08","Moffat_09","Moffat_10","Moffat_11","Moffat_15_S1")], genes=data[,1:2] );
g <- rep(c("MC38.OVA_Atg12.1_control","MC38.OVA_Atg12.1_IFNg"),each=3);
f <- factor( g );

keep.exprs <- filterByExpr( dge, group=f );
dge <- dge[keep.exprs,,keep.lib.sizes=FALSE];
dge <- calcNormFactors( dge,method="TMM" );

design <- model.matrix( ~ f );

y <- voom( dge, design );
fit <- lmFit( y, design );
fit2 <- eBayes( fit );

ttable <- topTable( fit2, coef="fMC38.OVA_Atg12.1_IFNg",number=Inf );
write.table( ttable,"topTable_Comp3_MC38.OVA_Atg12.1_IFNg_vs_Ctrl.txt",sep="\t",row.names=F );

dt <- decideTests( fit2 );
pdf("scatterplot_mc38Ova_Atg12.1_IFNg_Tx_18July19.pdf");
bigpar();
plotMD( fit2,column=2,status=dt[,2],main="MC38OVA - Atg12-1 sgRNA - +/- IFNg");
with(subset(ttable,adj.P.Val <= 0.05),text(AveExpr,logFC,labels=EntrezGene.ID,pos=4,cex=0.7));
abline(h=0,col="darkgrey");
dev.off();



#  4. Interaction - Intergenic +/- IFNγ vs. Atg12-1 +/-IFNγ
dge <- DGEList( counts=data[,c("Moffat_01","Moffat_02","Moffat_03","Moffat_04","Moffat_05","Moffat_06","Moffat_07","Moffat_08","Moffat_09","Moffat_10","Moffat_11","Moffat_15_S1")], genes=data[,1:2] );
sg <- factor(rep(c("Intergenic","Atg12.1"),each=6),levels=c("Intergenic","Atg12.1"));
tx <- factor(rep(c("Control","IFNg"),each=3,times=2));

keep.exprs <- filterByExpr( dge, group=paste(sg,tx,sep=".") );
dge <- dge[keep.exprs,,keep.lib.sizes=FALSE];
dge <- calcNormFactors( dge,method="TMM" );

design <- model.matrix( ~ sg * tx );

y <- voom( dge, design );
fit <- lmFit( y, design );
fit2 <- eBayes( fit );

ttable <- topTable( fit2, coef="sgAtg12.1:txIFNg",number=Inf );
write.table( ttable,"topTable_Comp4_MC38.OVA_Atg12.1_IFNg_vs_Ctrl_InteractionTerm.txt",sep="\t",row.names=F );

dt <- decideTests( fit2 );
pdf("scatterplot_mc38Ova_Atg12.1_IFNg_Interaction_18July19.pdf");
bigpar();
plotMD( fit2,column=4,status=dt[,2],main="MC38OVA - Atg12-1 sgRNA - +/- IFNg - Interaction");
with(subset(ttable,adj.P.Val <= 0.05),text(AveExpr,logFC,labels=EntrezGene.ID,pos=4,cex=0.7));
abline(h=0,col="darkgrey");
dev.off();





#  5. Intergenic control [14, 15, 16] vs. Atg12-1 control [20, 21, 22]
dge <- DGEList( counts=data[,c("Moffat_13","Moffat_14","Moffat_12","Moffat_19","Moffat_20","Moffat_21")], genes=data[,1:2] );
g <- rep(c("Renca-HA_Ctrl","Renca-HA_Atg12.1_control"),each=3);
f <- factor( g,levels=c("Renca-HA_Ctrl","Renca-HA_Atg12.1_control") );

keep.exprs <- filterByExpr( dge, group=f );
dge <- dge[keep.exprs,,keep.lib.sizes=FALSE];
dge <- calcNormFactors( dge,method="TMM" );

design <- model.matrix( ~ f );

y <- voom( dge, design );
fit <- lmFit( y, design );
fit2 <- eBayes( fit );

ttable <- topTable( fit2, coef="fRenca-HA_Atg12.1_control",number=Inf );
write.table( ttable,"topTable_Comp1_Renca-HA_Atg12.1_Ctrl.txt",sep="\t",row.names=F );

dt <- decideTests( fit2 );
pdf("scatterplot_rencaHa_Atg12.1_Ctrl_18July19.pdf");
bigpar();
plotMD( fit2,column=2,status=dt[,2],main="Renca-HA - Atg12.1 Control");
with(subset(ttable,adj.P.Val <= 0.05),text(AveExpr,logFC,labels=EntrezGene.ID,pos=4,cex=0.7));
abline(h=0,col="darkgrey");
dev.off();





#  6. Intergenic control #2 [38, 39, 40] vs. Fitm2-1 control [44, 45, 46]
dge <- DGEList( counts=data[,c("Moffat_37","Moffat_38","Moffat_39","Moffat_43","Moffat_44","Moffat_45")], genes=data[,1:2] );
g <- rep(c("RencaHA.Ctrl","Fitm2.Ctrl"),each=3);
f <- factor( g,levels=c("RencaHA.Ctrl","Fitm2.Ctrl") );

keep.exprs <- filterByExpr( dge, group=f );
dge <- dge[keep.exprs,,keep.lib.sizes=FALSE];
dge <- calcNormFactors( dge,method="TMM" );

design <- model.matrix( ~ f );

y <- voom( dge, design );
fit <- lmFit( y, design );
fit2 <- eBayes( fit );

ttable <- topTable( fit2, coef="fFitm2.Ctrl",number=Inf );
write.table( ttable,"topTable_Comp6_RencaHA_Fitm2_Effect.txt",sep="\t",row.names=F );

dt <- decideTests( fit2 );
pdf("scatterplot_RencaHA_Fitm2_Effect_18July19.pdf");
bigpar();
plotMD( fit2,column=2,status=dt[,2],main="Renca-HA - Fitm2 Effect");
with(subset(ttable,adj.P.Val <= 0.05),text(AveExpr,logFC,labels=EntrezGene.ID,pos=4,cex=0.7));
abline(h=0,col="darkgrey");
dev.off();



#  7. Intergenic control [14, 15, 16] vs. Intergenic + IFNγ [17, 18, 19]
dge <- DGEList( counts=data[,c("Moffat_13","Moffat_14","Moffat_12","Moffat_16","Moffat_17","Moffat_18")], genes=data[,1:2] );
g <- rep(c("RencaHA.Ctrl","RencaHA.IFNg"),each=3);
f <- factor( g );

keep.exprs <- filterByExpr( dge, group=f );
dge <- dge[keep.exprs,,keep.lib.sizes=FALSE];
dge <- calcNormFactors( dge,method="TMM" );

design <- model.matrix( ~ f );

y <- voom( dge, design );
fit <- lmFit( y, design );
fit2 <- eBayes( fit );

ttable <- topTable( fit2, coef="fRencaHA.IFNg",number=Inf );
write.table( ttable,"topTable_Comp7_RencaHA_IFNg_Effect",sep="\t",row.names=F );

dt <- decideTests( fit2 );
pdf("scatterplot_RencaHA_intergenic_IFNg_Tx_18July19.pdf");
bigpar();
plotMD( fit2,column=2,status=dt[,2],main="Renca-HA +/- IFNg");
with(subset(ttable,adj.P.Val <= 0.05),text(AveExpr,logFC,labels=EntrezGene.ID,pos=4,cex=0.7));
abline(h=0,col="darkgrey");
dev.off();



#  8. Atg12-1 control [20, 21, 22] vs. Atg12-1 IFNγ [23, 24, 25]
dge <- DGEList( counts=data[,c("Moffat_19","Moffat_20","Moffat_21","Moffat_22","Moffat_23","Moffat_24")], genes=data[,1:2] );
g <- rep(c("RencaHA.Atg12","RencaHA.Atg12.IFNg"),each=3);
f <- factor( g );

keep.exprs <- filterByExpr( dge, group=f );
dge <- dge[keep.exprs,,keep.lib.sizes=FALSE];
dge <- calcNormFactors( dge,method="TMM" );

design <- model.matrix( ~ f );

y <- voom( dge, design );
fit <- lmFit( y, design );
fit2 <- eBayes( fit );

ttable <- topTable( fit2, coef="fRencaHA.Atg12.IFNg",number=Inf );
write.table( ttable,"topTable_Comp8_RencaHA_Atg12.1_IFNg_Effect",sep="\t",row.names=F );

dt <- decideTests( fit2 );
pdf("scatterplot_RencaHA_Atg12_IFNg_Tx_18July19.pdf");
bigpar();
plotMD( fit2,column=2,status=dt[,2],main="Renca-HA.Atg12 +/- IFNg");
with(subset(ttable,adj.P.Val <= 0.05),text(AveExpr,logFC,labels=EntrezGene.ID,pos=4,cex=0.7));
abline(h=0,col="darkgrey");
dev.off();



#  9. Intergenic control #2 [38, 39, 40] vs. Intergenic + IFNγ #2 [41, 42, 43]
dge <- DGEList( counts=data[,c("Moffat_37","Moffat_38","Moffat_39","Moffat_40","Moffat_41","Moffat_42")], genes=data[,1:2] );
g <- rep(c("RencaHA.2.Ctrl","RencaHA.2.IFNg"),each=3);
f <- factor( g );

keep.exprs <- filterByExpr( dge, group=f );
dge <- dge[keep.exprs,,keep.lib.sizes=FALSE];
dge <- calcNormFactors( dge,method="TMM" );

design <- model.matrix( ~ f );

y <- voom( dge, design );
fit <- lmFit( y, design );
fit2 <- eBayes( fit );

ttable <- topTable( fit2, coef="fRencaHA.2.IFNg",number=Inf );
write.table( ttable,"topTable_Comp9_RencaHA.2_IFNg_Effect",sep="\t",row.names=F );

dt <- decideTests( fit2 );
pdf("scatterplot_RencaHA.2_intergenic_IFNg_Tx_18July19.pdf");
bigpar();
plotMD( fit2,column=2,status=dt[,2],main="Renca-HA#2 +/- IFNg");
with(subset(ttable,adj.P.Val <= 0.05),text(AveExpr,logFC,labels=EntrezGene.ID,pos=4,cex=0.7));
abline(h=0,col="darkgrey");
dev.off();



#  10. Fitm2-1 control [44, 45, 46] vs. Fitm2-1 + IFNγ [47, 48, 49]
dge <- DGEList( counts=data[,c("Moffat_43","Moffat_44","Moffat_45","Moffat_46","Moffat_47","Moffat_48")], genes=data[,1:2] );
g <- rep(c("RencaHA.Fitm2","RencaHA.Fitm2.IFNg"),each=3);
f <- factor( g );

keep.exprs <- filterByExpr( dge, group=f );
dge <- dge[keep.exprs,,keep.lib.sizes=FALSE];
dge <- calcNormFactors( dge,method="TMM" );

design <- model.matrix( ~ f );

y <- voom( dge, design );
fit <- lmFit( y, design );
fit2 <- eBayes( fit );

ttable <- topTable( fit2, coef="fRencaHA.Fitm2.IFNg",number=Inf );
write.table( ttable,"topTable_Comp10_RencaHA_Fitm2_IFNg_Effect",sep="\t",row.names=F );

dt <- decideTests( fit2 );
pdf("scatterplot_RencaHA_FItm2_IFNg_Tx_18July19.pdf");
bigpar();
plotMD( fit2,column=2,status=dt[,2],main="Renca-HA.Fitm2 +/- IFNg");
with(subset(ttable,adj.P.Val <= 0.05),text(AveExpr,logFC,labels=EntrezGene.ID,pos=4,cex=0.7));
abline(h=0,col="darkgrey");
dev.off();






#  11. RENCA-HA Intergenic +/- IFNγ vs. Atg12-1 +/-IFNγ
dge <- DGEList( counts=data[,c("Moffat_13","Moffat_14","Moffat_12","Moffat_16","Moffat_17","Moffat_18","Moffat_19","Moffat_20","Moffat_21","Moffat_22","Moffat_23","Moffat_24")], genes=data[,1:2] );
sg <- factor(rep(c("Intergenic","Atg12.1"),each=6),levels=c("Intergenic","Atg12.1"));
tx <- factor(rep(c("Control","IFNg"),each=3,times=2));

keep.exprs <- filterByExpr( dge, group=paste(sg,tx,sep=".") );
dge <- dge[keep.exprs,,keep.lib.sizes=FALSE];
dge <- calcNormFactors( dge,method="TMM" );

design <- model.matrix( ~ sg * tx );

y <- voom( dge, design );
fit <- lmFit( y, design );
fit2 <- eBayes( fit );

ttable.atg <- topTable( fit2, coef="sgAtg12.1",number=Inf );
ttable.tx <- topTable( fit2, coef="txIFNg",number=Inf );
ttable.interaction <- topTable( fit2, coef="sgAtg12.1:txIFNg",number=Inf );
write.table( ttable.atg,"topTable_RencaHA_Atg12.1_Effect.txt",sep="\t",row.names=F );
write.table( ttable.tx,"topTable_RencaHA_IFNg_Effect.txt",sep="\t",row.names=F );
write.table( ttable.interaction,"topTable_RencaHA_Atg12.1_IFNg_InteractionTerm.txt",sep="\t",row.names=F );

dt <- decideTests( fit2 );
pdf("scatterplot_RencaHA_Atg12.1_IFNg_Interaction_18July19.pdf");
bigpar();
plotMD( fit2,column=4,status=dt[,4],main="RencaHA - Atg12-1 sgRNA - +/- IFNg - Interaction");
with(subset(ttable.interaction,adj.P.Val <= 0.05),text(AveExpr,logFC,labels=EntrezGene.ID,pos=4,cex=0.7));
abline(h=0,col="darkgrey");
dev.off();

pdf("scatterplot_RencaHA_Atg12-1_sgRNA_Effect_18July19.pdf");
bigpar();
plotMD( fit2,column=2,status=dt[,2],main="RencaHA - Atg12-1 sgRNA Effect");
with(subset(ttable.atg,adj.P.Val <= 0.05),text(AveExpr,logFC,labels=EntrezGene.ID,pos=4,cex=0.7));
abline(h=0,col="darkgrey");
dev.off();

pdf("scatterplot_RencaHA_IFNg_Effect_18July19.pdf");
bigpar();
plotMD( fit2,column=3,status=dt[,3],main="RencaHA - IFNg Effect");
with(subset(ttable.tx,adj.P.Val <= 0.05),text(AveExpr,logFC,labels=EntrezGene.ID,pos=4,cex=0.7));
abline(h=0,col="darkgrey");
dev.off();




#  12. RENCA-HA Intergenic +/- IFNγ #2 vs. Fitm2-1 +/-IFNγ
dge <- DGEList( counts=data[,c("Moffat_37","Moffat_38","Moffat_39","Moffat_40","Moffat_41","Moffat_42","Moffat_43","Moffat_44","Moffat_45","Moffat_46","Moffat_47","Moffat_48")], genes=data[,1:2] );
sg <- factor(rep(c("Intergenic","Fitm2"),each=6),levels=c("Intergenic","Fitm2"));
tx <- factor(rep(c("Control","IFNg"),each=3,times=2));

keep.exprs <- filterByExpr( dge, group=paste(sg,tx,sep=".") );
dge <- dge[keep.exprs,,keep.lib.sizes=FALSE];
dge <- calcNormFactors( dge,method="TMM" );

design <- model.matrix( ~ sg * tx );

y <- voom( dge, design );
fit <- lmFit( y, design );
fit2 <- eBayes( fit );

ttable.fitm2 <- topTable( fit2, coef="sgFitm2",number=Inf );
ttable.tx <- topTable( fit2, coef="txIFNg",number=Inf );
ttable.interaction <- topTable( fit2, coef="sgFitm2:txIFNg",number=Inf );
write.table( ttable.fitm2,"topTable_RencaHA.2_Fitm2_Effect.txt",sep="\t",row.names=F );
write.table( ttable.tx,"topTable_RencaHA.2_IFNg_Effect.txt",sep="\t",row.names=F );
write.table( ttable.interaction,"topTable_RencaHA.2_Fitm2_IFNg_InteractionTerm.txt",sep="\t",row.names=F );

dt <- decideTests( fit2 );
pdf("scatterplot_RencaHA.2_Fitm2_IFNg_Interaction_18July19.pdf");
bigpar();
plotMD( fit2,column=4,status=dt[,4],main="RencaHA#2 - Fitm2 sgRNA - +/- IFNg - Interaction");
with(subset(ttable.interaction,adj.P.Val <= 0.05),text(AveExpr,logFC,labels=EntrezGene.ID,pos=4,cex=0.7));
abline(h=0,col="darkgrey");
dev.off();

pdf("scatterplot_RencaHA.2_Fitm2_sgRNA_Effect_18July19.pdf");
bigpar();
plotMD( fit2,column=2,status=dt[,2],main="RencaHA#2 - Fitm2 sgRNA Effect");
with(subset(ttable.atg,adj.P.Val <= 0.05),text(AveExpr,logFC,labels=EntrezGene.ID,pos=4,cex=0.7));
abline(h=0,col="darkgrey");
dev.off();

pdf("scatterplot_RencaHA.2_IFNg_Effect_18July19.pdf");
bigpar();
plotMD( fit2,column=3,status=dt[,3],main="RencaHA#2 - IFNg Effect");
with(subset(ttable.tx,adj.P.Val <= 0.05),text(AveExpr,logFC,labels=EntrezGene.ID,pos=4,cex=0.7));
abline(h=0,col="darkgrey");
dev.off();




#  13. Intergenic control [26, 27, 28] vs. Fitm2-1 control [32, 33, 34]
dge <- DGEList( counts=data[,c("Moffat_25","Moffat_26","Moffat_27","Moffat_31","Moffat_32","Moffat_33")], genes=data[,1:2] );
g <- rep(c("B16.Intergenic","B16.Fitm2.1"),each=3);
f <- factor( g,levels=c("B16.Intergenic","B16.Fitm2.1") );

keep.exprs <- filterByExpr( dge, group=f );
dge <- dge[keep.exprs,,keep.lib.sizes=FALSE];
dge <- calcNormFactors( dge,method="TMM" );

design <- model.matrix( ~ f );

y <- voom( dge, design );
fit <- lmFit( y, design );
fit2 <- eBayes( fit );

ttable <- topTable( fit2, coef="fB16.Fitm2.1",number=Inf );
write.table( ttable,"topTable_Comp13_B16_Fitm2.1_vs_Intergenic.txt",sep="\t",row.names=F );

dt <- decideTests( fit2 );
pdf("scatterplot_B16_Fitm2.1_vs_Intergenic_18July19.pdf");
bigpar();
plotMD( fit2,column=2,status=dt[,2],main="B16 Fitm2 vs Intergenic");
with(subset(ttable,adj.P.Val <= 0.05),text(AveExpr,logFC,labels=EntrezGene.ID,pos=4,cex=0.7));
abline(h=0,col="darkgrey");
dev.off();



#  14. Intergenic control [26, 27, 28] vs. Intergenic + IFNγ [29, 30, 31]
dge <- DGEList( counts=data[,c("Moffat_25","Moffat_26","Moffat_27","Moffat_28","Moffat_29","Moffat_30")], genes=data[,1:2] );
g <- rep(c("B16.Intergenic","B16.IFNg"),each=3);
f <- factor( g,levels=c("B16.Intergenic","B16.IFNg") );

keep.exprs <- filterByExpr( dge, group=f );
dge <- dge[keep.exprs,,keep.lib.sizes=FALSE];
dge <- calcNormFactors( dge,method="TMM" );

design <- model.matrix( ~ f );

y <- voom( dge, design );
fit <- lmFit( y, design );
fit2 <- eBayes( fit );

ttable <- topTable( fit2, coef="fB16.IFNg",number=Inf );
write.table( ttable,"topTable_Comp14_B16_Intergenic_IFNg.txt",sep="\t",row.names=F );

dt <- decideTests( fit2 );
pdf("scatterplot_B16_Intergenic_IFNg_18July19.pdf");
bigpar();
plotMD( fit2,column=2,status=dt[,2],main="B16 Intergenic + IFNg");
with(subset(ttable,adj.P.Val <= 0.05),text(AveExpr,logFC,labels=EntrezGene.ID,pos=4,cex=0.7));
abline(h=0,col="darkgrey");
dev.off();





#  15. Fitm2-1 control [32, 33, 34] vs. Fitm2-1 + IFNγ [35, 36, 37]
dge <- DGEList( counts=data[,c("Moffat_31","Moffat_32","Moffat_33","Moffat_34","Moffat_35","Moffat_36")], genes=data[,1:2] );
g <- rep(c("B16.Fitm2","B16.Fitm2.IFNg"),each=3);
f <- factor( g,levels=c("B16.Fitm2","B16.Fitm2.IFNg") );

keep.exprs <- filterByExpr( dge, group=f );
dge <- dge[keep.exprs,,keep.lib.sizes=FALSE];
dge <- calcNormFactors( dge,method="TMM" );

design <- model.matrix( ~ f );

y <- voom( dge, design );
fit <- lmFit( y, design );
fit2 <- eBayes( fit );

ttable <- topTable( fit2, coef="fB16.Fitm2.IFNg",number=Inf );
write.table( ttable,"topTable_Comp15_B16_Fitm2_IFNg.txt",sep="\t",row.names=F );

dt <- decideTests( fit2 );
pdf("scatterplot_B16_Fitm2_IFNg_18July19.pdf");
bigpar();
plotMD( fit2,column=2,status=dt[,2],main="B16 Fitm2 + IFNg");
with(subset(ttable,adj.P.Val <= 0.05),text(AveExpr,logFC,labels=EntrezGene.ID,pos=4,cex=0.7));
abline(h=0,col="darkgrey");
dev.off();




#  16. B16 Intergenic +/- IFNγ vs. Fitm2-1 +/-IFNγ
dge <- DGEList( counts=data[,paste("Moffat_",25:36,sep="")], genes=data[,1:2] );
sg <- factor(rep(c("Intergenic","Fitm2"),each=6),levels=c("Intergenic","Fitm2"));
tx <- factor(rep(c("Control","IFNg"),each=3,times=2));

keep.exprs <- filterByExpr( dge, group=paste(sg,tx,sep=".") );
dge <- dge[keep.exprs,,keep.lib.sizes=FALSE];
dge <- calcNormFactors( dge,method="TMM" );

design <- model.matrix( ~ sg * tx );

y <- voom( dge, design );
fit <- lmFit( y, design );
fit2 <- eBayes( fit );

ttable.fitm2 <- topTable( fit2, coef="sgFitm2",number=Inf );
ttable.tx <- topTable( fit2, coef="txIFNg",number=Inf );
ttable.interaction <- topTable( fit2, coef="sgFitm2:txIFNg",number=Inf );
write.table( ttable.fitm2,"topTable_B16_Fitm2_Effect.txt",sep="\t",row.names=F );
write.table( ttable.tx,"topTable_B16_IFNg_Effect.txt",sep="\t",row.names=F );
write.table( ttable.interaction,"topTable_B16_Fitm2_IFNg_InteractionTerm.txt",sep="\t",row.names=F );

dt <- decideTests( fit2 );
pdf("scatterplot_B16_Fitm2_IFNg_Interaction_18July19.pdf");
bigpar();
plotMD( fit2,column=4,status=dt[,4],main="B16 - Fitm2 sgRNA - +/- IFNg - Interaction");
with(subset(ttable.interaction,adj.P.Val <= 0.05),text(AveExpr,logFC,labels=EntrezGene.ID,pos=4,cex=0.7));
abline(h=0,col="darkgrey");
dev.off();

pdf("scatterplot_B16_Fitm2_sgRNA_Effect_18July19.pdf");
bigpar();
plotMD( fit2,column=2,status=dt[,2],main="B16 - Fitm2 sgRNA Effect");
with(subset(ttable.atg,adj.P.Val <= 0.05),text(AveExpr,logFC,labels=EntrezGene.ID,pos=4,cex=0.7));
abline(h=0,col="darkgrey");
dev.off();

pdf("scatterplot_B16_IFNg_Effect_18July19.pdf");
bigpar();
plotMD( fit2,column=3,status=dt[,3],main="B16 - IFNg Effect");
with(subset(ttable.tx,adj.P.Val <= 0.05),text(AveExpr,logFC,labels=EntrezGene.ID,pos=4,cex=0.7));
abline(h=0,col="darkgrey");
dev.off();





#  17. Renca-ATCC WT control 12hr [50, 51, 52] vs. Atg12 control 12hr [62, 63, 64]
dge <- DGEList( counts=data[,c("Moffat_49","Moffat_50","Moffat_51","Moffat_61","Moffat_62","Moffat_63")], genes=data[,1:2] );
g <- rep(c("Renca.ATCC.WT","Renca.ATCC.Atg12"),each=3);
f <- factor( g,levels=c("Renca.ATCC.WT","Renca.ATCC.Atg12") );

keep.exprs <- filterByExpr( dge, group=f );
dge <- dge[keep.exprs,,keep.lib.sizes=FALSE];
dge <- calcNormFactors( dge,method="TMM" );

design <- model.matrix( ~ f );

y <- voom( dge, design );
fit <- lmFit( y, design );
fit2 <- eBayes( fit );

ttable <- topTable( fit2, coef="fRenca.ATCC.Atg12",number=Inf );
write.table( ttable,"topTable_Comp17_Renca-ATCC_12hr_Atg12_vs_WT.txt",sep="\t",row.names=F );

dt <- decideTests( fit2 );
pdf("scatterplot_Renca-ATCC_12hr_Atg12_vs_WT_18July19.pdf");
bigpar();
plotMD( fit2,column=2,status=dt[,2],main="Renca-ATCC - 12h Atg12 vs WT");
with(subset(ttable,adj.P.Val <= 0.05),text(AveExpr,logFC,labels=EntrezGene.ID,pos=4,cex=0.7));
abline(h=0,col="darkgrey");
dev.off();



#  18. WT control 48hr [53, 54, 55] vs. Atg12 control 48hr [65, 66, 67]
dge <- DGEList( counts=data[,c("Moffat_52","Moffat_53","Moffat_54","Moffat_64","Moffat_65","Moffat_66")], genes=data[,1:2] );
g <- rep(c("Renca.ATCC.WT","Renca.ATCC.Atg12"),each=3);
f <- factor( g,levels=c("Renca.ATCC.WT","Renca.ATCC.Atg12") );

keep.exprs <- filterByExpr( dge, group=f );
dge <- dge[keep.exprs,,keep.lib.sizes=FALSE];
dge <- calcNormFactors( dge,method="TMM" );

design <- model.matrix( ~ f );

y <- voom( dge, design );
fit <- lmFit( y, design );
fit2 <- eBayes( fit );

ttable <- topTable( fit2, coef="fRenca.ATCC.Atg12",number=Inf );
write.table( ttable,"topTable_Comp18_Renca-ATCC_48h_Atg12_vs_WT.txt",sep="\t",row.names=F );

dt <- decideTests( fit2 );
pdf("scatterplot_Renca-ATCC_48h_Atg12_vs_WT_18July19.pdf");
bigpar();
plotMD( fit2,column=2,status=dt[,2],main="Renca-ATCC - 48h Atg12 vs WT");
with(subset(ttable,adj.P.Val <= 0.05),text(AveExpr,logFC,labels=EntrezGene.ID,pos=4,cex=0.7));
abline(h=0,col="darkgrey");
dev.off();




#  19. WT control 48hr [53, 54, 55] vs. Fitm2 control 48hr [74, 75, 76]
dge <- DGEList( counts=data[,c("Moffat_52","Moffat_53","Moffat_54","Moffat_73","Moffat_74","Moffat_75")], genes=data[,1:2] );
g <- rep(c("Renca.ATCC.WT","Renca.ATCC.Fitm2"),each=3);
f <- factor( g,levels=c("Renca.ATCC.WT","Renca.ATCC.Fitm2") );

keep.exprs <- filterByExpr( dge, group=f );
dge <- dge[keep.exprs,,keep.lib.sizes=FALSE];
dge <- calcNormFactors( dge,method="TMM" );

design <- model.matrix( ~ f );

y <- voom( dge, design );
fit <- lmFit( y, design );
fit2 <- eBayes( fit );

ttable <- topTable( fit2, coef="fRenca.ATCC.Fitm2",number=Inf );
write.table( ttable,"topTable_Comp19_Renca-ATCC_48h_Fitm2_vs_WT.txt",sep="\t",row.names=F );

dt <- decideTests( fit2 );
pdf("scatterplot_Renca-ATCC_48h_Fitm2_vs_WT_18July19.pdf");
bigpar();
plotMD( fit2,column=2,status=dt[,2],main="Renca-ATCC - 48h Fitm2 vs WT");
with(subset(ttable,adj.P.Val <= 0.05),text(AveExpr,logFC,labels=EntrezGene.ID,pos=4,cex=0.7));
abline(h=0,col="darkgrey");
dev.off();




#  20. WT control 12hr [50, 51, 52] vs. WT + TNF 12hr [59, 60, 61]
dge <- DGEList( counts=data[,c("Moffat_49","Moffat_50","Moffat_51","Moffat_58","Moffat_59","Moffat_60")], genes=data[,1:2] );
g <- rep(c("Renca.ATCC.WT","Renca.ATCC.Tnf"),each=3);
f <- factor( g,levels=c("Renca.ATCC.WT","Renca.ATCC.Tnf") );

keep.exprs <- filterByExpr( dge, group=f );
dge <- dge[keep.exprs,,keep.lib.sizes=FALSE];
dge <- calcNormFactors( dge,method="TMM" );

design <- model.matrix( ~ f );

y <- voom( dge, design );
fit <- lmFit( y, design );
fit2 <- eBayes( fit );

ttable <- topTable( fit2, coef="fRenca.ATCC.Tnf",number=Inf );
write.table( ttable,"topTable_Comp20_Renca-ATCC_12h_TNFtx_vs_WT.txt",sep="\t",row.names=F );

dt <- decideTests( fit2 );
pdf("scatterplot_Renca-ATCC_12h_TNFtx_vs_WT_18July19.pdf");
bigpar();
plotMD( fit2,column=2,status=dt[,2],main="Renca-ATCC - 12h +/- TNF");
with(subset(ttable,adj.P.Val <= 0.05),text(AveExpr,logFC,labels=EntrezGene.ID,pos=4,cex=0.7));
abline(h=0,col="darkgrey");
dev.off();




#  21. Atg12 control 12hr [62, 63, 64] vs. Atg12 + TNF 12hr [71, 72, 73]
dge <- DGEList( counts=data[,c("Moffat_61","Moffat_62","Moffat_63","Moffat_70","Moffat_71","Moffat_72")], genes=data[,1:2] );
g <- rep(c("Renca.ATCC.Atg12","Renca.ATCC.Atg12.Tnf"),each=3);
f <- factor( g,levels=c("Renca.ATCC.Atg12","Renca.ATCC.Atg12.Tnf") );

keep.exprs <- filterByExpr( dge, group=f );
dge <- dge[keep.exprs,,keep.lib.sizes=FALSE];
dge <- calcNormFactors( dge,method="TMM" );

design <- model.matrix( ~ f );

y <- voom( dge, design );
fit <- lmFit( y, design );
fit2 <- eBayes( fit );

ttable <- topTable( fit2, coef="fRenca.ATCC.Atg12.Tnf",number=Inf );
write.table( ttable,"topTable_Comp21_Renca-ATCC_12h_Atg12_TNFtx.txt",sep="\t",row.names=F );

dt <- decideTests( fit2 );
pdf("scatterplot_Renca-ATCC_12h_Atg12_TNFtx_18July19.pdf");
bigpar();
plotMD( fit2,column=2,status=dt[,2],main="Renca-ATCC - 12h Atg12 +/- TNF");
with(subset(ttable,adj.P.Val <= 0.05),text(AveExpr,logFC,labels=EntrezGene.ID,pos=4,cex=0.7));
abline(h=0,col="darkgrey");
dev.off();



#  22. WT control 48hr [53, 54, 55] vs. WT + IFN 48hr [56, 57, 58]
dge <- DGEList( counts=data[,c("Moffat_52","Moffat_53","Moffat_54","Moffat_55","Moffat_56","Moffat_57")], genes=data[,1:2] );
g <- rep(c("Renca.ATCC.WT","Renca.ATCC.WT.Ifn"),each=3);
f <- factor( g,levels=c("Renca.ATCC.WT","Renca.ATCC.WT.Ifn") );

keep.exprs <- filterByExpr( dge, group=f );
dge <- dge[keep.exprs,,keep.lib.sizes=FALSE];
dge <- calcNormFactors( dge,method="TMM" );

design <- model.matrix( ~ f );

y <- voom( dge, design );
fit <- lmFit( y, design );
fit2 <- eBayes( fit );

ttable <- topTable( fit2, coef="fRenca.ATCC.WT.Ifn",number=Inf );
write.table( ttable,"topTable_Comp22_Renca-ATCC_48h_WT_IFNtx.txt",sep="\t",row.names=F );

dt <- decideTests( fit2 );
pdf("scatterplot_Renca-ATCC_48h_WT_IFNtx_18July19.pdf");
bigpar();
plotMD( fit2,column=2,status=dt[,2],main="Renca-ATCC - 48h WT +/- IFNg");
with(subset(ttable,adj.P.Val <= 0.05),text(AveExpr,logFC,labels=EntrezGene.ID,pos=4,cex=0.7));
abline(h=0,col="darkgrey");
dev.off();



#  23. Atg12 control 48hr [65, 66, 67] vs. Atg12 + IFN 48hr [68, 69, 70]
dge <- DGEList( counts=data[,c("Moffat_64","Moffat_65","Moffat_66","Moffat_67","Moffat_68","Moffat_69")], genes=data[,1:2] );
g <- rep(c("Renca.ATCC.Atg","Renca.ATCC.Atg.Ifn"),each=3);
f <- factor( g,levels=c("Renca.ATCC.Atg","Renca.ATCC.Atg.Ifn") );

keep.exprs <- filterByExpr( dge, group=f );
dge <- dge[keep.exprs,,keep.lib.sizes=FALSE];
dge <- calcNormFactors( dge,method="TMM" );

design <- model.matrix( ~ f );

y <- voom( dge, design );
fit <- lmFit( y, design );
fit2 <- eBayes( fit );

ttable <- topTable( fit2, coef="fRenca.ATCC.Atg.Ifn",number=Inf );
write.table( ttable,"topTable_Comp23_Renca-ATCC_48h_Atg12_IFNtx.txt",sep="\t",row.names=F );

dt <- decideTests( fit2 );
pdf("scatterplot_Renca-ATCC_48h_Atg_IFNtx_18July19.pdf");
bigpar();
plotMD( fit2,column=2,status=dt[,2],main="Renca-ATCC - 48h Atg12 +/- IFNg");
with(subset(ttable,adj.P.Val <= 0.05),text(AveExpr,logFC,labels=EntrezGene.ID,pos=4,cex=0.7));
abline(h=0,col="darkgrey");
dev.off();



#  24. Fitm2-1 control 48hr [74, 75, 76] vs. Fitm2-1 + IFN 48hr[77, 78, 79]
dge <- DGEList( counts=data[,c("Moffat_73","Moffat_74","Moffat_75","Moffat_76","Moffat_77","Moffat_78")], genes=data[,1:2] );
g <- rep(c("Renca.ATCC.Fitm2","Renca.ATCC.Fitm2.Ifn"),each=3);
f <- factor( g,levels=c("Renca.ATCC.Fitm2","Renca.ATCC.Fitm2.Ifn") );

keep.exprs <- filterByExpr( dge, group=f );
dge <- dge[keep.exprs,,keep.lib.sizes=FALSE];
dge <- calcNormFactors( dge,method="TMM" );

design <- model.matrix( ~ f );

y <- voom( dge, design );
fit <- lmFit( y, design );
fit2 <- eBayes( fit );

ttable <- topTable( fit2, coef="fRenca.ATCC.Fitm2.Ifn",number=Inf );
write.table( ttable,"topTable_Comp24_Renca-ATCC_48h_Fitm2_IFNtx.txt",sep="\t",row.names=F );

dt <- decideTests( fit2 );
pdf("scatterplot_Renca-ATCC_48h_Fitm2_IFNtx_18July19.pdf");
bigpar();
plotMD( fit2,column=2,status=dt[,2],main="Renca-ATCC - 48h Fitm2 +/- IFNg");
with(subset(ttable,adj.P.Val <= 0.05),text(AveExpr,logFC,labels=EntrezGene.ID,pos=4,cex=0.7));
abline(h=0,col="darkgrey");
dev.off();



#  25. Renca-ATCC WT +/- TNF 12hr vs. Atg12 +/- TNF 12hr
dge <- DGEList( counts=data[,paste("Moffat_",c(49,50,51,58,59,60,61,62,63,70,71,72),sep="")], genes=data[,1:2] );
sg <- factor(rep(c("Intergenic","Atg12"),each=6),levels=c("Intergenic","Atg12"));
tx <- factor(rep(c("Control","TNF"),each=3,times=2));

keep.exprs <- filterByExpr( dge, group=paste(sg,tx,sep=".") );
dge <- dge[keep.exprs,,keep.lib.sizes=FALSE];
dge <- calcNormFactors( dge,method="TMM" );

design <- model.matrix( ~ sg * tx );

y <- voom( dge, design );
fit <- lmFit( y, design );
fit2 <- eBayes( fit );

ttable.atg <- topTable( fit2, coef="sgAtg12",number=Inf );
ttable.tx <- topTable( fit2, coef="txTNF",number=Inf );
ttable.interaction <- topTable( fit2, coef="sgAtg12:txTNF",number=Inf );
write.table( ttable.atg,"topTable_RencaATCC_Atg12_Effect.txt",sep="\t",row.names=F );
write.table( ttable.tx,"topTable_RencaATCC_12hTNF_Effect.txt",sep="\t",row.names=F );
write.table( ttable.interaction,"topTable_RencaATCC_Atg12_TNF_InteractionTerm.txt",sep="\t",row.names=F );

dt <- decideTests( fit2 );
pdf("scatterplot_RencaATCC_Atg12_TNF_Interaction_18July19.pdf");
bigpar();
plotMD( fit2,column=4,status=dt[,4],main="Renca-ATCC - Atg12 sgRNA - +/- TNF - Interaction");
with(subset(ttable.interaction,adj.P.Val <= 0.05),text(AveExpr,logFC,labels=EntrezGene.ID,pos=4,cex=0.7));
abline(h=0,col="darkgrey");
dev.off();

pdf("scatterplot_RencaATCC_Atg12_sgRNA_Effect_18July19.pdf");
bigpar();
plotMD( fit2,column=2,status=dt[,2],main="Renca-ATCC - Atg12 sgRNA Effect");
with(subset(ttable.atg,adj.P.Val <= 0.05),text(AveExpr,logFC,labels=EntrezGene.ID,pos=4,cex=0.7));
abline(h=0,col="darkgrey");
dev.off();

pdf("scatterplot_RencaATCC_TNF_Effect_18July19.pdf");
bigpar();
plotMD( fit2,column=3,status=dt[,3],main="Renca-ATCC - TNF Effect");
with(subset(ttable.tx,adj.P.Val <= 0.05),text(AveExpr,logFC,labels=EntrezGene.ID,pos=4,cex=0.7));
abline(h=0,col="darkgrey");
dev.off();




#  26. RencaATCC WT +/- IFN 48hr vs. Atg12 +/- IFN 48hr
dge <- DGEList( counts=data[,paste("Moffat_",c(52:57,64:69),sep="")], genes=data[,1:2] );
sg <- factor(rep(c("Intergenic","Atg12"),each=6),levels=c("Intergenic","Atg12"));
tx <- factor(rep(c("Control","IFNg"),each=3,times=2));

keep.exprs <- filterByExpr( dge, group=paste(sg,tx,sep=".") );
dge <- dge[keep.exprs,,keep.lib.sizes=FALSE];
dge <- calcNormFactors( dge,method="TMM" );

design <- model.matrix( ~ sg * tx );

y <- voom( dge, design );
fit <- lmFit( y, design );
fit2 <- eBayes( fit );

ttable.atg12 <- topTable( fit2, coef="sgAtg12",number=Inf );
ttable.tx <- topTable( fit2, coef="txIFNg",number=Inf );
ttable.interaction <- topTable( fit2, coef="sgAtg12:txIFNg",number=Inf );
write.table( ttable.atg12,"topTable_RencaATCC_48h_Atg12_Effect.txt",sep="\t",row.names=F );
write.table( ttable.tx,"topTable_RencaATCC_48h_IFNg_Effect.txt",sep="\t",row.names=F );
write.table( ttable.interaction,"topTable_RencaATCC_48h_Atg12_IFNg_InteractionTerm.txt",sep="\t",row.names=F );

dt <- decideTests( fit2 );
pdf("scatterplot_RencaATCC_48h_Atg12_IFNg_Interaction_18July19.pdf");
bigpar();
plotMD( fit2,column=4,status=dt[,4],main="Renca-ATCC - 48h Atg12 sgRNA - +/- IFNg - Interaction");
with(subset(ttable.interaction,adj.P.Val <= 0.05),text(AveExpr,logFC,labels=EntrezGene.ID,pos=4,cex=0.7));
abline(h=0,col="darkgrey");
dev.off();

pdf("scatterplot_RencaATCC_48h_Atg12_sgRNA_Effect_18July19.pdf");
bigpar();
plotMD( fit2,column=2,status=dt[,2],main="Renca-ATCC - 48h Atg12 sgRNA Effect");
with(subset(ttable.atg,adj.P.Val <= 0.05),text(AveExpr,logFC,labels=EntrezGene.ID,pos=4,cex=0.7));
abline(h=0,col="darkgrey");
dev.off();

pdf("scatterplot_RencaATCC_48h_IFNg_Effect_18July19.pdf");
bigpar();
plotMD( fit2,column=3,status=dt[,3],main="Renca-ATCC - 48h IFNg Effect");
with(subset(ttable.tx,adj.P.Val <= 0.05),text(AveExpr,logFC,labels=EntrezGene.ID,pos=4,cex=0.7));
abline(h=0,col="darkgrey");
dev.off();



#  27. WT +/- IFN 48hr vs. Fitm2 +/- IFN 48hr
dge <- DGEList( counts=data[,paste("Moffat_",c(52:57,73:78),sep="")], genes=data[,1:2] );
sg <- factor(rep(c("Intergenic","Fitm2"),each=6),levels=c("Intergenic","Fitm2"));
tx <- factor(rep(c("Control","IFNg"),each=3,times=2));

keep.exprs <- filterByExpr( dge, group=paste(sg,tx,sep=".") );
dge <- dge[keep.exprs,,keep.lib.sizes=FALSE];
dge <- calcNormFactors( dge,method="TMM" );

design <- model.matrix( ~ sg * tx );

y <- voom( dge, design );
fit <- lmFit( y, design );
fit2 <- eBayes( fit );

ttable.fitm2 <- topTable( fit2, coef="sgFitm2",number=Inf );
ttable.tx <- topTable( fit2, coef="txIFNg",number=Inf );
ttable.interaction <- topTable( fit2, coef="sgFitm2:txIFNg",number=Inf );
write.table( ttable.fitm2,"topTable_RencaATCC_48h_Fitm2_Effect.txt",sep="\t",row.names=F );
write.table( ttable.tx,"topTable_RencaATCC_48h_IFNg_Effect.txt",sep="\t",row.names=F );
write.table( ttable.interaction,"topTable_RencaATCC_48h_Fitm2_IFNg_InteractionTerm.txt",sep="\t",row.names=F );

dt <- decideTests( fit2 );
pdf("scatterplot_RencaATCC_48h_Fitm2_IFNg_Interaction_18July19.pdf");
bigpar();
plotMD( fit2,column=4,status=dt[,4],main="Renca-ATCC - 48h Fitm2 sgRNA - +/- IFNg - Interaction");
with(subset(ttable.interaction,adj.P.Val <= 0.05),text(AveExpr,logFC,labels=EntrezGene.ID,pos=4,cex=0.7));
abline(h=0,col="darkgrey");
dev.off();

pdf("scatterplot_RencaATCC_48h_Fitm2_sgRNA_Effect_18July19.pdf");
bigpar();
plotMD( fit2,column=2,status=dt[,2],main="Renca-ATCC - 48h Fitm2 sgRNA Effect");
with(subset(ttable.atg,adj.P.Val <= 0.05),text(AveExpr,logFC,labels=EntrezGene.ID,pos=4,cex=0.7));
abline(h=0,col="darkgrey");
dev.off();