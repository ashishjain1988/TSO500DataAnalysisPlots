#rm(list=ls())

#setwd("~/Documents/saksena_oncoplot")

library(openxlsx)
library(maftools)
library(trackViewer)

### Read in data
xl_file="~/Downloads/Landscape of somatic alterations in FCL-LFGT..relevant genes for oncoplot (1).xlsx"
oncoplot_data <- openxlsx::read.xlsx(xl_file, sheet="fake_maf")

oncoplot_data$Variant_Classification[oncoplot_data$Variant_Classification == "(splicing;UTR5)"]<-"Splicing_UTR5"
oncoplot_data$Variant_Classification[oncoplot_data$Variant_Classification == "(splicing)"]<-"Splicing"
oncoplot_data$Variant_Classification[oncoplot_data$Variant_Classification == "FS substitution"]<-"FS_Substitution"
oncoplot_data$Variant_Classification[oncoplot_data$Variant_Classification == "nonFS substitution"]<-"Non_FS_Substitution"
oncoplot_data$Variant_Classification[oncoplot_data$Variant_Classification == "stopgain"]<-"Stopgain"

variant_type_colors <- c("SNV"="#377EB8",
                         "FS_Substitution"="#4DAF4A",
                         "Non_FS_Substitution"="#ff008c",
                         "Stopgain"="#A65628",
                         "Splicing"="#FF7F00","Splicing_UTR5"="#ad7aff")
variant_types <- names(variant_type_colors)

mymaf <- read.maf(oncoplot_data, vc_nonSyn=variant_types)
mymaf@data$Hugo_Symbol
library(berryFunctions)

allgenes <- unique(mymaf@data[,c("Hugo_Symbol","Transcript")])
allgenes$Transcript <- unlist(lapply(strsplit(allgenes$Transcript,"\\."), "[[", 1))
out_dir="~/Downloads/"
if(!dir.exists(out_dir)){dir.create(out_dir, recursive = T)}
pdf(file.path(out_dir,"lollipop_plots-PFAM.pdf"),width = 10,height = 5)
genes<-c("TNFRSF14","B2M","XIAP","CARD11","PTPRD","SOCS1")
for (mygene in genes) {
  # mygene="TNFRSF14"
  try(
    lollipopPlot(mymaf, mygene,showMutationRate=F,labelPos = NULL,showLegend = T, domainAlpha=1,showDomainLabel=F, refSeqID = allgenes$Transcript[allgenes$Hugo_Symbol==mygene],AACol="AAChange",colors=variant_type_colors, roundedRect = T),
    silent = T

  )
}
dev.off()


