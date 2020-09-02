# Date: Sep 28, 2018
# Author: Shirley Hui, Ruth Isserlin
# Takes supplied list of genes and uses gProfiler to perform enrichment analysis to determine which pathways summarize the genes supplied.
# Input: list of core killing genes (Fig 1G)
# Output: Pathway enrichment output file
if(!"RCy3" %in% installed.packages()){
  install.packages("BiocManager")
  BiocManager::install("RCy3")
}
library(RCy3)
library(RCurl)

tryCatch(expr = { library("gProfileR")}, 
         error = function(e) { install.packages("gProfileR")}, finally = library("gProfileR"))

# Function to run gprofiler using the gprofiler library
# The function takes the returned gprofiler results and formats it to the generic EM input file
# function returns a data frame in the generic EM file format.
runGprofiler <- function(genes,current_organism = "mmusculus", ###mmusculus for mouse!!
                         significant_only = T, set_size_max = 200, 
                         set_size_min = 3, filter_gs_size_min = 5 , exclude_iea = F){
  
  gprofiler_results <- gprofiler(genes ,
                                 significant=significant_only,ordered_query = F,
                                 exclude_iea=exclude_iea,max_set_size = set_size_max,
                                 min_set_size = set_size_min,
                                 correction_method = "fdr",
                                 organism = current_organism,
                                 src_filter = c("GO:BP","GO:MF","GO:CC","REAC","KEGG","WP","TF","MIRNA","HPA","CORUM"))
  
  # Filter results
  gprofiler_results <- gprofiler_results[which(gprofiler_results[,'term.size'] >= 3
                                               & gprofiler_results[,'overlap.size'] >= filter_gs_size_min ),]
  
  # gProfileR returns corrected p-values only.  Set p-value to corrected p-value
  if(dim(gprofiler_results)[1] > 0){
    
    gprofiler_results_filename <-"/Users/Keith/Desktop/Nature Final Revision/Revision Docs/Pathway analysis/Output/gprofiler_results/coreCTLgenes_gprofiler.txt"
    write.table(gprofiler_results,gprofiler_results_filename,col.name=TRUE,sep="\t",row.names=FALSE,quote=FALSE)
    
    em_results <- cbind(gprofiler_results[,
                                          c("term.id","term.name","p.value","p.value")], 1,
                        gprofiler_results[,"intersection"])
    colnames(em_results) <- c("Name","Description", "pvalue","qvalue","phenotype","genes")
    
    return(em_results)
  } else {
    return("no gprofiler results for supplied query")
  }
}

  genes = read.delim("/Users/Keith/Desktop/Nature Final Revision/Revision Docs/Pathway analysis/Input/coreCTLgenes.txt", header = F)
genes = as.vector(genes[,1])
gprofiler_results = runGprofiler(genes)

# Write out the g:Profiler results
em_results_filename <-"/Users/Keith/Desktop/Nature Final Revision/Revision Docs/Pathway analysis/Output/em_file/coreCTLgenes_sim0.7.txt"
write.table(gprofiler_results,em_results_filename,col.name=TRUE,sep="\t",row.names=FALSE,quote=FALSE)

em_command = paste('enrichmentmap build analysisType="generic" ', 
                   'pvalue=',"0.001", 'qvalue=',"0.05",
                   'similaritycutoff=',"0.7",
                   'coeffecients=',"JACCARD+OVERLAP",
                   'enrichmentsDataset1=',em_results_filename ,
                   sep=" ")

# Enrichment map command will return the suid of newly created network.
em_network_suid <- commandsRun(em_command)
renameNetwork("Cluster1_enrichmentmap", network=as.numeric(em_network_suid))
