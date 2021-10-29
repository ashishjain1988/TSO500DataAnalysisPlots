library(openxlsx)
library(circlize)
library(ComplexHeatmap)
library(RColorBrewer)
library(data.table)
library(fs)
library(stringr)
library(maftools)

#SNV data format
xl_file="~/Downloads/Landscape of somatic alterations in FCL-LFGT..relevant genes for oncoplot (1).xlsx"
oncoplot_data <- openxlsx::read.xlsx(xl_file, sheet="fake_maf")
snvs<-openxlsx::read.xlsx("~/Downloads/Landscape of somatic alterations in FCL-LFGT..relevant genes for oncoplot.xlsx",sheet = 1)
patientCaseMapping<-unique(snvs[,c(2,4)])
piMap<-patientCaseMapping$Case.number
names(piMap)<-patientCaseMapping$Patient.ID
oncoplot_data$Variant_Classification[oncoplot_data$Variant_Classification == "(splicing;UTR5)"]<-"Splicing UTR5"
oncoplot_data$Variant_Classification[oncoplot_data$Variant_Classification == "(splicing)"]<-"Splicing"
oncoplot_data$Variant_Classification[oncoplot_data$Variant_Classification == "FS substitution"]<-"FS Substitution"
oncoplot_data$Variant_Classification[oncoplot_data$Variant_Classification == "nonFS substitution"]<-"Non FS Substitution"
oncoplot_data$Variant_Classification[oncoplot_data$Variant_Classification == "stopgain"]<-"Stopgain"
mutation_colors <- c("Splicing_UTR5"="#ad7aff",Missense_Mutation="black",SNV="#377EB8","FS_Substitution"="#4DAF4A","Non_FS_Substitution"="#ff008c",Splicing="#FF7F00",Stopgain="#A65628",Multi_Hit="#FFFF33",Duplication="#5abad8",Del="darkred",no_variants="#d6d6d6", Pathogenic="black",VUS="grey50")
names(mutation_colors) <- gsub("_"," ",names(mutation_colors))
variant_types <- names(mutation_colors)
oncoplot_data$Tumor_Sample_Barcode<-piMap[oncoplot_data$Tumor_Sample_Barcode]
#mymaf <- read.maf(oncoplot_data, vc_nonSyn=variant_types)

###CNV data format
cnvs<-read.table("~/Downloads/CNV data for case6-5152.txt",sep=",",header = TRUE)
custom.cnv.data = data.frame(Gene = cnvs$Gene,Sample_name = cnvs$Patient.ID,CN = cnvs$Dup.del,stringsAsFactors = FALSE)
custom.cnv.data$CN[custom.cnv.data$CN == "DUP"]<-"Duplication"
custom.cnv.data$Sample_name <-piMap[custom.cnv.data$Sample_name]

##Load clinical features
f<-openxlsx::read.xlsx("~/Downloads/Landscape of somatic alterations in FCL-LFGT..relevant genes for oncoplot.xlsx", sheet=2)[1:9,]
clinicalFeatures<-data.frame(Tumor_Sample_Barcode=f$Case.number,IG_Gene_Rearr=trimws(f$IG),BCL2=f$`BCL2.t(14;18)`)
maf.filtered = read.maf(maf = oncoplot_data,cnTable = custom.cnv.data,clinicalData = clinicalFeatures,verbose = FALSE,vc_nonSyn=variant_types)

##Colors for clinical annotations
clin_data<-maf.filtered@clinical.data[,c("Tumor_Sample_Barcode","IG_Gene_Rearr","BCL2")]
IG_Gene_RearrColor <- c("#78DD9E","#A54EE1","#D2B4C8")
names(IG_Gene_RearrColor)<-c("Clonal","Indeterminate","Polyclonal")
BCL2_color<-c("brown")
names(BCL2_color)<-c("Neg")
clin_data_colors <- list(IG_Gene_Rearr=IG_Gene_RearrColor,BCL2=BCL2_color)

##Loading pathway information
pathwayAnnot<-openxlsx::read.xlsx("~/Downloads/PathwayInformation.xlsx",sheet = 1)
pathwayAnnot$Gene<-trimws(pathwayAnnot$Gene)
pathwayAnnot1<-merge(data.frame(Gene=maf.filtered@data$Hugo_Symbol),pathwayAnnot,by="Gene",all.x=TRUE)
pathwayAnnot1$Pathway[is.na(pathwayAnnot1$Pathway)]<-"Other"


###Code to plot oncoplot for TSO500 data
require(ComplexHeatmap)
### Structure info about the fraction of the cohort that has each gene mutated
frac_mut <- data.frame(Hugo_Symbol=maf.filtered@gene.summary$Hugo_Symbol,
                       frac_mut=(maf.filtered@gene.summary$AlteredSamples/as.numeric(maf.filtered@summary$summary[3])),
                       #mutation_count=maf.filtered@gene.summary$total + maf.filtered@gene.summary$CNV_total,
                       mutation_count=maf.filtered@gene.summary$total,
                       #mutation_count=maf.filtered@gene.summary$CNV_total,
                       stringsAsFactors = F)

frac_mut <- frac_mut %>% dplyr::filter(!(Hugo_Symbol %in% c("HLA-A","HLA-B","HLA-C")))
#frac_mut <- frac_mut %>% dplyr::filter(!(Hugo_Symbol %in% c("HLA-A","HLA-B","HLA-C")) & (mutation_count > 1))

ngene_max=100
target_frac = sort(frac_mut$frac_mut, decreasing = T)[min(ngene_max,nrow(frac_mut))]
cohort_freq_thresh = 0.01
cohort_freq_thresh <- max(c(cohort_freq_thresh,target_frac))
### Select genes based on the frequency threshold
#frac_mut <- frac_mut[order(frac_mut$frac_mut,frac_mut$mutation_count,decreasing = T),]
frac_mut <- frac_mut[order(frac_mut$frac_mut,decreasing = T),]
freq_genes <- frac_mut$Hugo_Symbol[frac_mut$frac_mut >= cohort_freq_thresh]
freq_genes <- freq_genes[1:min(ngene_max,length(freq_genes))]
if (length(freq_genes) == 0) {
  stop("No genes to plot; change the frequency threshold to include more genes.")
}
if (length(freq_genes) > 100) {
  target_frac = round(sort(frac_mut$frac_mut, decreasing = T)[min(ngene_max,nrow(frac_mut))],2)
  # stop(paste0("Too many genes for oncoplot. Trying setting the cohort mutated fraction to > ", target_frac))
  warning(paste0("Too many genes for oncoplot. Trying setting the cohort mutated fraction to > ", target_frac))
  # return(NA)
}
gene_list <- list(freq_genes)
reasons <- paste0("Cohort Freq > ",round(cohort_freq_thresh,digits = 3))

### Collect genes to plot
genes_for_oncoplot <- data.frame(Hugo_Symbol=c(), reason=c())
for (i in 1:length(gene_list)) {
  if (is.na(gene_list[[i]][1])) {
    next
  }
  genes_for_oncoplot <- rbind(genes_for_oncoplot,
                              data.frame(Hugo_Symbol=gene_list[[i]],
                                         reason=reasons[i]))
}
genes_for_oncoplot <- cbind(genes_for_oncoplot,
                            frac=frac_mut$frac_mut[match(genes_for_oncoplot$Hugo_Symbol, frac_mut$Hugo_Symbol)],
                            reason1=pathwayAnnot1$Pathway[match(genes_for_oncoplot$Hugo_Symbol, pathwayAnnot1$Gene)])

genes_for_oncoplot <- genes_for_oncoplot[order(genes_for_oncoplot$reason1, -genes_for_oncoplot$frac),]
# browser()
### Split the oncoplot based on the reason for picking the gene
###   Here, we're only picked based on the frequency
###   But this framework is useful for plotting genes picked using various criteria
split_idx=genes_for_oncoplot$reason1
split_colors <- rainbow(length(unique(split_idx)))
# names(split_colors) <- as.character(genes_for_oncoplot$reason[!duplicated(genes_for_oncoplot$reason)])
names(split_colors) <- unique(split_idx)
split_colors <- list(Reason=split_colors)
split_idx <- as.factor(split_idx)
split_idx<-factor(split_idx,levels=c("Immune modulation","Chromatin remodeling","Apoptosis","Cell Growth/differentiation",
                                     "JAK-STAT pathway","Ras-MAPK/Ras-PI3K pathways","Other signaling pathways","Cell cycle/cytoskeleton","Transcription factor","B-cell development","Other"))

# source("scripts/helper_functions.oncoplot.R")
### Make matrix to plot, and order it correctly
#print(genes_for_oncoplot$Hugo_Symbol)
oncomat <- createOncoMatrix(maf.filtered, g=genes_for_oncoplot$Hugo_Symbol, add_missing = F)$oncoMatrix
oncomat <- oncomat[match(genes_for_oncoplot$Hugo_Symbol,rownames(oncomat)), ]
onco_genes <- rownames(oncomat)
oncomat.plot <- oncomat

### Set the height of the plot based on number of genes
onco_height=NULL
if (is.null(onco_height)) {
  onco_height=max(round(0.2*nrow(oncomat.plot),0),5)
}

### Make the mutation type names prettier by removing the underscore
# my_mut_col <- mutation_colors
# names(mutation_colors) <- gsub("_"," ",names(mutation_colors))
oncomat.plot <- gsub("_"," ",oncomat.plot)

myanno=NULL
if (!is.null(clin_data)) {
  # browser()
  anno_data <- data.frame(clin_data[match(colnames(oncomat.plot), clin_data$Tumor_Sample_Barcode),],stringsAsFactors = F)
  row.names(anno_data) <- anno_data$Tumor_Sample_Barcode
  anno_data <- anno_data[,!colnames(anno_data) %in% "Tumor_Sample_Barcode", drop=F]
  if (ncol(anno_data) > 0) {
    ###Make changes in the text for the annotation
    colnames(anno_data)<-c("IG Gene\nRearrangement","BCL2 rearrangement\nby FISH")
    myanno <- HeatmapAnnotation(df=anno_data,col = list("IG Gene\nRearrangement"=clin_data_colors$IG_Gene_Rearr,"BCL2 rearrangement\nby FISH"=clin_data_colors$BCL2),
                                simple_anno_size = unit(10, "mm"),
                                annotation_name_gp =  gpar(fontsize = 12,fontface = 2),
                                annotation_legend_param=list(labels_gp = gpar(fontsize = 14),title_gp = gpar(fontsize = 16, fontface = 2),
                                                             nrow = 4,
                                                             legend_direction = "vertical"))
  }
}

## Show total burden for top annotation
variant_type_data <- data.frame(maf.filtered@variant.classification.summary)
rownames(variant_type_data) <- variant_type_data$Tumor_Sample_Barcode
colnames(variant_type_data) <- gsub("_"," ",colnames(variant_type_data))
#variant_type_data <- variant_type_data[,c(-1,-ncol(variant_type_data))]
variant_type_data <- variant_type_data[,!names(variant_type_data)%in%c("Tumor Sample Barcode","total","CNV total")]
variant_type_data <- variant_type_data[match(colnames(oncomat.plot), rownames(variant_type_data)),
                                       rev(order(colSums(variant_type_data)))]
# browser()
var_anno_colors <- mutation_colors[match(colnames(variant_type_data), names(mutation_colors))]
###Make changes in the text/data for the top mutation histogram
# top_ha = HeatmapAnnotation("# of Mutations" = anno_barplot(variant_type_data, gp = gpar(fill = var_anno_colors), border = F,rot=90,height = unit(4, "cm")),
#                            #"Samples" = anno_text(sort(colnames(oncomat.plot)),rot=0),
#                            annotation_name_side = "left",annotation_name_gp = gpar(fontsize = 12),gap=unit(1,"mm"))
top_ha = HeatmapAnnotation("TMB (Mut/MB)" = anno_barplot(f$TMB[as.numeric(rownames(variant_type_data))], gp = gpar(fill = "mediumorchid4"), border = F,rot=90,height = unit(4, "cm")),
                           #"Samples" = anno_text(sort(colnames(oncomat.plot)),rot=0),
                           annotation_name_side = "left",annotation_name_gp = gpar(fontsize = 12),gap=unit(1,"mm"))

# browser()

pct_anno <- paste0(prettyNum(frac_mut$frac_mut[match(onco_genes, frac_mut$Hugo_Symbol)]*100,digits=1),"%")
# OG_TSG_DF1<-data.frame(Gene=onco_genes)
# OG_TSG_DF <- merge(OG_TSG_DF1,pathdbDF,by="Gene",all.x=TRUE)
# OG_TSG_DF$Pathway[is.na(OG_TSG_DF$Pathway)]<-"NA"
# OG_TSG_DF$OG_TSG[is.na(OG_TSG_DF$OG_TSG)]<-"NA"
# oncogene_anno <- paste0(OG_TSG_DF$OG_TSG[match(onco_genes, OG_TSG_DF$Gene)])
#left_ha = rowAnnotation("Oncogene"=anno_text(oncogene_anno,gp = gpar(fontsize = 12)),"Cohort Pct"=anno_text(pct_anno,gp = gpar(fontsize = 12)), show_annotation_name=F,gap=unit(1,"mm"))
right_ha = rowAnnotation("% Mutated Sample"=anno_text(pct_anno,gp = gpar(fontsize = 14)), show_annotation_name=TRUE,gap=unit(1,"mm"))
#row.names(oncomat.plot)<-paste0(row.names(oncomat.plot)," (",pct_anno,")")
left_ha = rowAnnotation("Pathways"=anno_text(genes_for_oncoplot$reason1,gp = gpar(fontsize = 12),just = "left"), show_annotation_name=TRUE)
#pathwayAnnotation <- rowAnnotation(month = anno_text(month.name[1:10], just = "center",location = unit(0.5, "npc"), show_name = TRUE),annotation_name_rot = 0)

# print(oncomat.plot)
### Make the oncoplot
# alter_fun$Deletion<-function(x, y, w, h) {
#   grid.rect(x, y, w-unit(0.25, "mm"), h*0.33,
#             gp = gpar(fill = mutation_colors["Del"], col = NA))
# }
# alter_fun$Amplification<-function(x, y, w, h) {
#   grid.rect(x, y, w-unit(0.25, "mm"), h*0.33,
#             gp = gpar(fill = mutation_colors["Amp"], col = NA))
# }
col = mutation_colors
onco_base_default <- oncoPrint(oncomat.plot,get_type = function(x) strsplit(x, ";")[[1]], alter_fun = function(x, y, w, h, v) {
  #print(v)
  n = sum(v)  # how many alterations for current gene in current sample
  h = h*0.9
  # use `names(which(v))` to correctly map between `v` and `col`
  if(n) {grid.rect(x, y - h*0.5 + 1:n/n*h, w*0.95, 1/n*h,
                  gp = gpar(fill = col[names(which(v))], col = NA), just = "top")}
  else
  {
      grid.rect(x, y, w-unit(0.5, "mm"), h,
      gp = gpar(fill = "#CCCCCC", col = NA))
  }
},
                               col=mutation_colors,
                               row_order=1:nrow(oncomat.plot),
                               name="oncoplot",
                               column_order = sort(colnames(oncomat.plot)),
                               show_pct = F,
                               right_annotation = right_ha,
                               row_split=split_idx,
                               row_title = NULL,
                               bottom_annotation = myanno,
                               top_annotation = top_ha,
                               left_annotation = left_ha,
                               show_column_names = TRUE,
                               column_names_gp = gpar(fontsize = 20,fontface = 2),
                               column_names_rot = 0,
                               #alter_fun_is_vectorized = T),
                               row_names_gp = gpar(fontsize = 14),
                               heatmap_legend_param = list(title = "Alterations",title_gp = gpar(fontsize = 16, fontface = 2),labels_gp = gpar(fontsize = 14)))#,

pdf(paste0("~/Downloads/TSO500-Oncoplot-test.pdf"),width = 14,height = 12)
#ComplexHeatmap::draw(g, show_annotation_legend = TRUE)
draw(onco_base_default, show_annotation_legend = TRUE)
dev.off()

# dataForCbioportal<-data.frame(Sample_ID=piMap[snvs$Patient.ID],Cancer_Type=rep("PCFCL",nrow(snvs)),Chromosome=snvs$Chr,Start_Position=snvs$Start,End_Position=snvs$End,Reference_Allele=snvs$Ref,Variant_Allele=snvs$Alt)
# write.table(dataForCbioportal,"~/Downloads/landscape.cbioportal.tsv",sep = "\t",row.names = FALSE,col.names = TRUE,quote = F)

# custom.snv.data <- data.frame(Gene = snvs$Gene,Sample_name = snvs$Patient.ID,CN = snvs$Exonic.function,stringsAsFactors = FALSE)
# custom.snv.data$CN[custom.snv.data$CN == "(splicing;UTR5)"]<-"Splicing_UTR5"
# custom.snv.data$CN[custom.snv.data$CN == "(splicing)"]<-"Splicing"
# custom.snv.data$CN[custom.snv.data$CN == "FS substitution"]<-"FS_Substitution"
# custom.snv.data$CN[custom.snv.data$CN == "nonFS substitution"]<-"Non_FS_Substitution"
# custom.snv.data$CN[custom.snv.data$CN == "stopgain"]<-"Stopgain"

# custom.cn.data <- custom.cnv.data
# custom.cn.data$Sample_name <-piMap[custom.cn.data$Sample_name]


createOncoMatrix = function(m, g = NULL, chatty = TRUE, add_missing = FALSE){

  if(is.null(g)){
    stop("Please provde atleast two genes!")
  }

  subMaf = subsetMaf(maf = m, genes = g, includeSyn = FALSE, mafObj = FALSE)

  if(nrow(subMaf) == 0){
    if(add_missing){
      numericMatrix = matrix(data = 0, nrow = length(g), ncol = length(levels(getSampleSummary(x = m)[,Tumor_Sample_Barcode])))
      rownames(numericMatrix) = g
      colnames(numericMatrix) = levels(getSampleSummary(x = m)[,Tumor_Sample_Barcode])

      oncoMatrix = matrix(data = "", nrow = length(g), ncol = length(levels(getSampleSummary(x = m)[,Tumor_Sample_Barcode])))
      rownames(oncoMatrix) = g
      colnames(oncoMatrix) = levels(getSampleSummary(x = m)[,Tumor_Sample_Barcode])

      vc = c("")
      names(vc) = 0

      return(list(oncoMatrix = oncoMatrix, numericMatrix = numericMatrix, vc = vc))
    }else{
      return(NULL)
    }
  }

  if(add_missing){
    subMaf[, Hugo_Symbol := factor(x = Hugo_Symbol, levels = g)]
  }

  cnv_events = c(c("Amp", "Del"), as.character(subMaf[Variant_Type == "CNV"][, .N, Variant_Classification][, Variant_Classification]))
  cnv_events = unique(cnv_events)

  oncomat = data.table::dcast(data = subMaf[,.(Hugo_Symbol, Variant_Classification, Tumor_Sample_Barcode)], formula = Hugo_Symbol ~ Tumor_Sample_Barcode,
                              fun.aggregate = function(x, cnv = cnv_events){
                                #x = unique(as.character(x)) #>=2 distinct variant classification = Multi_Hit
                                x = as.character(x) # >= 2 same/distinct variant classification = Multi_Hit See #347
                                xad = x[x %in% cnv]
                                xvc = x[!x %in% cnv]

                                if(length(xvc)>0){
                                  xvc = ifelse(test = length(xvc) > 1, yes = paste(xvc,collapse = ";"), no = xvc)
                                }

                                x = ifelse(test = length(xad) > 0, yes = paste(xad, xvc, sep = ';'), no = xvc)
                                x = gsub(pattern = ';$', replacement = '', x = x)
                                x = gsub(pattern = '^;', replacement = '', x = x)
                                return(x)
                              } , value.var = 'Variant_Classification', fill = '', drop = FALSE)

  #convert to matrix
  data.table::setDF(oncomat)
  rownames(oncomat) = oncomat$Hugo_Symbol
  oncomat = as.matrix(oncomat[,-1, drop = FALSE])

  variant.classes = as.character(unique(subMaf[,Variant_Classification]))
  variant.classes = c('',variant.classes, 'Multi_Hit')
  names(variant.classes) = 0:(length(variant.classes)-1)

  #Complex variant classes will be assigned a single integer.
  vc.onc = unique(unlist(apply(oncomat, 2, unique)))
  vc.onc = vc.onc[!vc.onc %in% names(variant.classes)]
  names(vc.onc) = rep(as.character(as.numeric(names(variant.classes)[length(variant.classes)])+1), length(vc.onc))
  variant.classes2 = c(variant.classes, vc.onc)

  oncomat.copy <- oncomat
  #Make a numeric coded matrix
  for(i in 1:length(variant.classes2)){
    oncomat[oncomat == variant.classes2[i]] = names(variant.classes2)[i]
  }

  #If maf has only one gene
  if(nrow(oncomat) == 1){
    mdf  = t(matrix(as.numeric(oncomat)))
    rownames(mdf) = rownames(oncomat)
    colnames(mdf) = colnames(oncomat)
    return(list(oncoMatrix = oncomat.copy, numericMatrix = mdf, vc = variant.classes))
  }

  #convert from character to numeric
  mdf = as.matrix(apply(oncomat, 2, function(x) as.numeric(as.character(x))))
  rownames(mdf) = rownames(oncomat.copy)


  #If MAF file contains a single sample, simple sorting is enuf.
  if(ncol(mdf) == 1){
    sampleId = colnames(mdf)
    mdf = as.matrix(mdf[order(mdf, decreasing = TRUE),])
    colnames(mdf) = sampleId

    oncomat.copy = as.matrix(oncomat.copy[rownames(mdf),])
    colnames(oncomat.copy) = sampleId

    return(list(oncoMatrix = oncomat.copy, numericMatrix = mdf, vc = variant.classes))
  } else{
    #Sort by rows as well columns if >1 samples present in MAF
    #Add total variants per gene
    mdf = cbind(mdf, variants = apply(mdf, 1, function(x) {
      length(x[x != "0"])
    }))
    #Sort by total variants
    mdf = mdf[order(mdf[, ncol(mdf)], decreasing = TRUE), ]
    #colnames(mdf) = gsub(pattern = "^X", replacement = "", colnames(mdf))
    nMut = mdf[, ncol(mdf)]

    mdf = mdf[, -ncol(mdf)]

    mdf.temp.copy = mdf #temp copy of original unsorted numeric coded matrix

    mdf[mdf != 0] = 1 #replacing all non-zero integers with 1 improves sorting (& grouping)
    tmdf = t(mdf) #transposematrix
    mdf = t(tmdf[do.call(order, c(as.list(as.data.frame(tmdf)), decreasing = TRUE)), ]) #sort

    mdf.temp.copy = mdf.temp.copy[rownames(mdf),] #organise original matrix into sorted matrix
    mdf.temp.copy = mdf.temp.copy[,colnames(mdf)]
    mdf = mdf.temp.copy

    #organise original character matrix into sorted matrix
    oncomat.copy <- oncomat.copy[,colnames(mdf)]
    oncomat.copy <- oncomat.copy[rownames(mdf),]

    return(list(oncoMatrix = oncomat.copy, numericMatrix = mdf, vc = variant.classes, cnvc = cnv_events))
  }
}

### List defining functions for color and shape of cells in oncoplot
# alter_fun = list(
#   background = function(x, y, w, h) {
#     grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
#               gp = gpar(fill = "#CCCCCC", col = NA))
#   },
#   # "0" = function(x, y, w, h) {
#   #   grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
#   #             gp = gpar(fill = "#CCCCCC", col = NA))
#   # },
#   "Nonsense Mutation" = function(x, y, w, h) {
#     grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
#               gp = gpar(fill = mutation_colors["Nonsense Mutation"], col = NA))
#   },
#   "Missense Mutation" = function(x, y, w, h) {
#     grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
#               gp = gpar(fill = mutation_colors["Missense Mutation"], col = NA))
#   },
#   "Frame Shift Del" = function(x, y, w, h) {
#     grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
#               gp = gpar(fill = mutation_colors["Frame Shift Del"], col = NA))
#   },
#   "In Frame Ins" = function(x, y, w, h) {
#     grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
#               gp = gpar(fill = mutation_colors["In Frame Ins"], col = NA))
#   },
#   "Splice Site" = function(x, y, w, h) {
#     grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
#               gp = gpar(fill = mutation_colors["Splice Site"], col = NA))
#   },
#   "Multi Hit" = function(x, y, w, h) {
#     grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
#               gp = gpar(fill = mutation_colors["Multi Hit"], col = NA))
#   },
#   "Frame Shift Ins" = function(x, y, w, h) {
#     grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
#               gp = gpar(fill = mutation_colors["Frame Shift Ins"], col = NA))
#   },
#   "In Frame Del" = function(x, y, w, h) {
#     grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
#               gp = gpar(fill = mutation_colors["In Frame Del"], col = NA))
#   },
#   "Nonstop Mutation" = function(x, y, w, h) {
#     grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
#               gp = gpar(fill = mutation_colors["Nonstop Mutation"], col = NA))
#   },
#   "Translation Start Site" = function(x, y, w, h) {
#     grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
#               gp = gpar(fill = mutation_colors["Translation Start Site"], col = NA))
#   },
#   "Duplication" = function(x, y, w, h) {
#     grid.rect(x, y, w-unit(0.25, "mm"), h*0.33,
#               gp = gpar(fill = mutation_colors["Duplication"], col = NA))
#   },
#   "Del" = function(x, y, w, h) {
#     grid.rect(x, y, w-unit(0.25, "mm"), h*0.33,
#               gp = gpar(fill = mutation_colors["Del"], col = NA))
#   },
#   "no variants" = function(x, y, w, h) {
#     grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
#               # gp = gpar(fill = "#e0e0e0", col = NA))
#               gp = gpar(fill = "#CCCCCC", col = NA))
#   },
#   "Pathogenic" = function(x, y, w, h) {
#     # grid.points(x, y, pch = 18, size=w, gp=gpar(col=col["pathogenic"]))
#     # grid.rect(x, y, w*0.7, h*0.2,
#     #           gp = gpar(fill = col["pathogenic"], col = NA))
#     # grid.rect(x, y, w*0.1, h*0.7,
#     #           gp = gpar(fill = col["pathogenic"], col = NA))
#     grid.rect(x, y, w*0.8, h*0.8,
#               gp = gpar(col = mutation_colors["Pathogenic"], fill = NA, lwd=5))
#   },
#   "VUS" = function(x, y, w, h) {
#     # grid.points(x, y, pch = 3, size=w,gp=gpar(col=col["VUS"], lwd=3))
#     # grid.rect(x, y, w*0.2, h-unit(0.5, "mm"),
#     #           gp = gpar(fill = col["VUS"], col = NA))
#     grid.rect(x, y, w*0.8, h*0.8,
#               gp = gpar(col = mutation_colors["VUS"], fill = NA, lwd=5))
#   }
#   ,
#   "Splicing UTR5" = function(x, y, w, h) {
#     grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
#               gp = gpar(fill = mutation_colors["Splicing UTR5"], col = NA))
#   }
#   ,
#   "FS Substitution" = function(x, y, w, h) {
#     grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
#               gp = gpar(fill = mutation_colors["FS Substitution"], col = NA))
#   }
#   ,
#   "Non FS Substitution" = function(x, y, w, h) {
#     grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
#               gp = gpar(fill = mutation_colors["Non FS Substitution"], col = NA))
#   }
#   ,
#   "Splicing" = function(x, y, w, h) {
#     grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
#               gp = gpar(fill = mutation_colors["Splicing"], col = NA))
#   }
#   ,
#   "SNV" = function(x, y, w, h) {
#     grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
#               gp = gpar(fill = mutation_colors["SNV"], col = NA))
#   }
#   ,
#   "Stopgain" = function(x, y, w, h) {
#     grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
#               gp = gpar(fill = mutation_colors["Stopgain"], col = NA))
#   }
# )
