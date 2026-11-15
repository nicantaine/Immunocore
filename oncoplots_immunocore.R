library(ComplexHeatmap)
library(tidyverse)

library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)

get_type_fun = function(x) strsplit(x, " ")[[1]]

col = c("DEL"="slategray2","Frame_Shift_Del"="plum2","Missense_Mutation"="red",
        "Frame_Shift_Ins" = "purple1","Splice_Site" = "navy","Nonsense_Mutation" = "black",
        "Nonstop_Mutation"="purple2","Translation_Start_Site" = "purple","Splice_Region"="blue","In_Frame_Del"="maroon3","In_Frame_Ins"="maroon4")

#removed for confidentiality 
pid<-c("")

#removed for confidentiality 
pid_by_tsb<-c("")

#removed for confidentiality 
ha_all = HeatmapAnnotation(Patient=pid,TIC = pid_by_tsb,
                           col = list(TIC=c(""),
                                      Patient=c("")), annotation_height = 5)

heatmap_legend_param = list(title = "Legend", at = c("Frame_Shift_Del","Frame_Shift_Ins","Missense_Mutation",
                                                     "Nonsense_Mutation","Splice_Site","Frame_Shift_Del","Translation_Start_Site",
                                                     "In_Frame_Del","In_Frame_Ins","Nonstop_Mutation"))

alter_fun <- list( background = function(x, y, w, h) {
  grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
            gp = gpar(fill = "#CCCCCC", col = NA))
},
Frame_Shift_Del = function(x, y, w, h) {
  grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
            gp = gpar(fill = col["Frame_Shift_Del"], col = NA))
},
Missense_Mutation = function(x, y, w, h) {
  grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
            gp = gpar(fill = col["Missense_Mutation"], col = NA))
},
Frame_Shift_Ins = function(x, y, w, h) {
  grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
            gp = gpar(fill = col["Frame_Shift_Ins"], col = NA))
},
Splice_Site = function(x, y, w, h) {
  grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
            gp = gpar(fill = col["Splice_Site"], col = NA))
},
Nonsense_Mutation = function(x, y, w, h) {
  grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
            gp = gpar(fill = col["Nonsense_Mutation"], col = NA))
},
Translation_Start_Site = function(x, y, w, h) {
  grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
            gp = gpar(fill = col["Translation_Start_Site"], col = NA))
},
Splice_Region = function(x, y, w, h) {
  grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
            gp = gpar(fill = col["Splice_Region"], col = NA))
},
Targeted_Region = function(x, y, w, h) {
  grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
            gp = gpar(fill = col["Targeted_Region"], col = NA))
},
In_Frame_Ins = function(x, y, w, h) {
  grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
            gp = gpar(fill = col["In_Frame_Ins"], col = NA))
},
In_Frame_Del = function(x, y, w, h) {
  grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
            gp = gpar(fill = col["In_Frame_Del"], col = NA))
},
Nonstop_Mutation = function(x, y, w, h) {
  grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
            gp = gpar(fill = col["Nonstop_Mutation"], col = NA))
})



cancer<-read.csv("", sep="\t",header=F)
column_title="All DNASeq Data"
colnames(cancer)<-c("Hugo_Symbol","Variant_Classification","Tumor_Sample_Barcode")
cancer <-cancer %>% group_by(Hugo_Symbol) %>% mutate(row=row_number()) %>% tidyr::pivot_wider(names_from = Tumor_Sample_Barcode,values_from = Variant_Classification)
cancer=cancer[-2]
cancer<-as.matrix(cancer)
cancer[is.na(cancer)] = ""
cancer<-as.data.frame(cancer)
cancer <- cancer %>% group_by(Hugo_Symbol) %>% summarise_all(funs(trimws(paste(., collapse = ' '))))
cancer<-as.data.frame(cancer)
rownames(cancer)<-cancer[,1]
drop<-"Hugo_Symbol"
cancer<-cancer[ , !(names(cancer) %in% drop)]
cancer<-as.matrix(cancer)

top_index <- order(apply(cancer, 1, function(x) sum(x != "")), decreasing = TRUE)[1:50]
mydata<-cancer[top_index,]
mydata<-mydata[,tsb_by_bor]

oncoPrint(mydata,alter_fun = alter_fun, col = col, column_title = "Top 50 Mutated Genes", 
          heatmap_legend_param = heatmap_legend_param, get_type = get_type_fun, 
          remove_empty_columns = TRUE,remove_empty_rows = TRUE,
          show_column_names = FALSE, row_names_gp = gpar(fontsize=6),
          bottom_annotation = ha_all, top_annotation = NULL)

