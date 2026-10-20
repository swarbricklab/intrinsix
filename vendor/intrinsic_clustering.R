#Intrinsic Cluster Data and Make Heatmap
#5/14/25
#Alec Lobanov
#Perou Lab
suppressPackageStartupMessages({
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(cluster)
  library(ctc)
  library(yaml)
})
source('modules/PAM50/clustering_scripts/clustering_scripts.R')
set.seed(1)

medianCtr<-function(x){
  annAll <- dimnames(x)
  medians <- apply(x,1,median,na.rm=T)
  x <- t(scale(t(x),center=medians,scale=F))
  dimnames(x) <- annAll
  return(x)
}

config <- yaml::read_yaml("config/PAM50/config.yaml")

#### Load in Metadata/add PAM50 ####
pam50_meta <- read.csv(config$intermediates$pam50_metadata,row.names = 1)
pam50scores <- read.table(paste0(config$outs$output_dir,'/PAM50scores_pam50scores.txt'),header = T,sep = '\t',row.names = 1,check.names = F)

#### Add PAM50 Scores to Heatmap Metadata ####
pam50_meta <- merge(pam50_meta,pam50scores["Call"],by = 0)
colnames(pam50_meta) <- c("Row.names","clinical_subtype","tissue_id","PAM50_Call")
rownames(pam50_meta) <- pam50_meta$Row.names
pam50_meta[,1] <- NULL

#### Identify Claudin Low Samples ####
intrinsic_expr_medCtr <- as.matrix(read.table(config$intermediates$int_expr_medCtr,header = T,row.names = 1,sep = '\t',check.names = F))

##Scripts from clustering scripts
cluster_files <- centroid_clustering(intrinsic_expr_medCtr,path = config$outs$int_heatmap_dir, write_files = T)
clow_samples <- identify_claudin_low(cluster_files,intrinsic_expr_medCtr)
pam50_meta$PAM50_Call_CLOW <- ifelse(rownames(pam50_meta) %in% clow_samples,'CLOW',pam50_meta$PAM50_Call)


#### Find Silhouette Width Scores ####
intrinsic_expr <- read.table(config$intermediates$int_expr,header = T,sep = ',',check.names = F,row.names = 1)

#Calculate Silhouette Width
clusters <- as.numeric(as.factor(pam50_meta$PAM50_Call_CLOW))
names(clusters) <- rownames(pam50_meta)
dist_matrix <- as.dist(1 - cor(intrinsic_expr, method = "pearson"))
silhouette_res <- silhouette(clusters,dist_matrix)
rownames(silhouette_res) <- rownames(pam50_meta)

## Sil Width Figure (Optional)
colors <- c('red','gold','pink','blue','lightblue','green')
pdf(paste0(config$outs$int_heatmap_dir,"/silhouette_plot.pdf"), width = 9, height = 8)
plot(silhouette_res,border = NA,col = colors)
abline(v=0.00)
dev.off()

## Add To metadata
colnames(silhouette_res) <- c('silhouette_cluster','neighbor','silhouette_width')
pam50_meta$silhouette_width <- silhouette_res[,3]
pam50_meta$sil_width_bin_0.00 <- ifelse(pam50_meta$silhouette_width < 0.00,"Negative","Positive")

#save final metadata
write.csv(pam50_meta,config$intermediates$pam50_metadata) #Overwrites it, can change if needed


#### Create Heatmap ####

#match metadata and data
if (!all.equal(rownames(pam50_meta),colnames(intrinsic_expr_medCtr))){
  pam50_meta <- pam50_meta[match(colnames(intrinsic_expr_medCtr), intrinsic_expr_medCtr),] #order info/gene list file to match matrix
}

#Heatmap Code
column_ha <- HeatmapAnnotation(
  "Clinical_subtype" = pam50_meta$clinical_subtype,
  "PAM50_zhao4+CLOW" = pam50_meta$PAM50_Call_CLOW,
  "Silhouette_width" = pam50_meta$sil_width_bin_0.00,
  col = list(
    'Clinical_subtype' = c("ERpos_HER2neg" = "lightblue",
                           "HER2pos_ERneg"="pink",
                           "HER2pos_ERpos" = "plum",
                           "ERneg_HER2neg" ="darkorange",
                           "TNBC" = "orange"),
    'PAM50_zhao4+CLOW' = c("Basal" = "red",
                           "CLOW" = "yellow",
                           "Her2" = "pink",
                           "LumA" = "darkblue",
                           "LumB" = "lightblue",
                           "Normal" = "green"),
    'Silhouette_width' = c("Negative"="#b2e2a8","Positive"='#fdbf6f')
  ))

#for reading in files from Cluster3.0
## Would love to get this straight from list outputs but need to think how to save it
my_cluster_columns <- xcluster2r(file = paste0(config$outs$int_heatmap_dir,list.files(path = config$outs$int_heatmap_dir,pattern = "\\.atr$")),
                                 distance = "pearson", labels = FALSE, fast = FALSE, clean = TRUE, dec='.')
my_cluster_rows <- xcluster2r(file = paste0(config$outs$int_heatmap_dir,list.files(path = config$outs$int_heatmap_dir,pattern = "\\.gtr$")),
                              distance = "pearson", labels = FALSE, fast = FALSE, clean = TRUE, dec='.')

#draw heatmap; lots of aesthetic options here
hm <- Heatmap(intrinsic_expr_medCtr, name = "Gene exp.",
              show_column_names = F,
              show_row_names = F,
              top_annotation = column_ha,
              row_title = "genes",
              column_title = "samples", column_names_side = "top",
              row_names_gp = gpar(fontsize = 8, fontface = "bold"),
              column_names_gp = gpar(fontsize = 8, fontface = "bold"),
              column_names_rot = 90, column_title_side = "bottom",
              cluster_rows = my_cluster_rows,
              cluster_columns = my_cluster_columns,
              col = colorRamp2(breaks = c(-4,0,4), c("blue","black", "yellow")), #sets the colors for your heatmap
              column_dend_height = unit(2, "cm"), #for extending the length of the array dendrogram
              row_dend_reorder = TRUE,
              column_dend_reorder = TRUE,
              column_names_max_height = unit(10, "cm"),
              use_raster = F #for forcing plot to be made in higher resolution
)


# Save heatmap
pdf(paste0(config$outs$int_heatmap_dir,"/intrinsic_medctr_heatmap.pdf"), width = 11, height = 8)
draw(hm, column_title = paste0("BRCA Major Atlas WTS Untreated Dataset,n=",nrow(pam50_meta)))
dev.off()