library(pheatmap)
library(tidyverse)
library(readxl)
library(rlist)
library(ggplot2)
require("ggplot2")
library(ggpubr)
library(ggthemes)
library(cowplot)

setwd('/research_jude/rgs01_jude/groups/xmagrp/projects/SingleCellDNA/common/DS-AML/WGS/tracking/Zhikai/figures/figure2_heatmap/')
input = read_excel("../DS_AML_zl.xlsx",sheet = "markers")
dna = input %>% select(SJID, Marker_hg38, nM_E1_WGS, nT_E1_WGS, nM_H1_WGS, nT_H1_WGS, nM_S1_WGS, nT_S1_WGS, Region, Class, AAChange, Gene)
genecount = dna %>% select(SJID, Gene, Class) %>% group_by(Gene) %>% summarise(n=n()) %>% arrange(desc(n))
genecount$Gene = factor(genecount$Gene, levels=genecount$Gene)
dna_input = dna %>% select(SJID, Gene, Class) %>% 
  mutate(newClass = ifelse(is.na(Class),"ExonDeletion",
                           ifelse(grepl("missense", Class), "missense",
                                  ifelse(grepl("proteinIns", Class), "proteinIns",
                                         ifelse(grepl("silent", Class), "silent",
                                                ifelse(grepl("frameshift", Class), "frameshift",
                                                       ifelse(grepl("startloss", Class), "startloss", Class)))))))
gata1 = dna_input %>% filter(Gene == "GATA1") %>% 
  group_by(newClass) %>% summarise(n=n()) %>%
  arrange(desc(n))
x = dna_input %>% arrange(match(newClass, gata1$newClass)) %>%
  filter(Gene == "GATA1")
myorder = c(as.character(gata1$newClass),"promoter")
dna_input$SJID = factor(dna_input$SJID, ordered=TRUE, 
       levels = unique(x$SJID))
dna_input$Gene = factor(dna_input$Gene, levels=rev(genecount$Gene))
gg = ggplot(dna_input, aes(x=SJID, y=Gene, fill=newClass, height=.8, width=.8)) + 
  geom_tile(color="white", position = "identity", linewidth=0.5)
gg = gg + theme_minimal()
gg = gg + labs(x=NULL, y=NULL, title="GATA1")
gg = gg + theme_tufte(base_family="Helvetica")
gg = gg + theme(axis.ticks=element_blank())+
  theme(axis.text.y = element_text(size = 12),
        axis.text.x=element_text(size=4))
gg = gg + theme(axis.text.x=element_text(angle = 90, hjust = 1)) +
  theme(legend.position="top")
gg
#dna_input$Gene = factor(dna_input$Gene, levels=genecount$Gene)
#dna_input$newClass = factor(dna_input$newClass, levels=c(myorder,"promoter"))
#rrange(factor(Reg, levels = LETTERS[c(3, 1, 2)])
###########
###########
###########
mylen = c()
for (gene in unique(genecount$Gene)){
  mylen = c(mylen, nchar(gene))
}
maxlen = max(mylen)
dna_input$newClass = factor(dna_input$newClass, 
                            levels = unique(dna_input$newClass))
figlst = list()
genelst = unique(genecount$Gene)
#genelst = c("GATA1","JAK3 ")
mycolors = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5")
myclass = as.character(unique(dna_input$newClass))
coltable = data.frame(cbind(mycolors, myclass))
for (gene in unique(genecount$Gene)){
  custom_breaks = 0.5
  thiscount = genecount %>% filter(Gene == gene) %>% pull(n)
  thisylabel = paste(gene," (",thiscount,")",sep="")
  custom_labels = thisylabel
  x = dna_input %>% filter(Gene == gene) %>% select(-Class) %>% 
    complete(SJID) %>% group_by(SJID, newClass) %>% summarise(n=n()) %>% 
    mutate(prop=ifelse(is.na(newClass),0,n/sum(n)))
  owncolor = coltable %>% filter(myclass %in% as.character(unique(x$newClass))) %>% pull(mycolors)
  #x$newClass = factor(x$newClass,levels = unique(dna_input$newClass))
  if (gene == "GATA1"){
    fig = ggplot(x, aes(fill=newClass, y=prop, x=SJID)) + 
      geom_bar(position="stack", stat="identity") +
      scale_y_continuous(breaks=custom_breaks, labels=custom_labels)+
      theme(axis.ticks.x = element_blank(),
            #axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_blank(),
            axis.title.x=element_blank(),legend.position = "top",
            plot.margin=unit(c(0,0,0,0),"cm")) +
      scale_fill_manual(values = owncolor)
  }else{
    fig = ggplot(x, aes(fill=newClass, y=prop, x=SJID)) + 
      geom_bar(position="stack", stat="identity") +
      scale_y_continuous(breaks=custom_breaks, labels=custom_labels)+
      theme(axis.ticks.x = element_blank(),
            #axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_blank(),
            axis.title.x=element_blank(),legend.position = "none",
            plot.margin=unit(c(0,0,0,0),"cm")) +
      scale_fill_manual(values = owncolor)
  }
  figlst = list.append(figlst, fig)
}
library(egg)
ggarrange(plots = figlst, nrow = length(unique(genecount$Gene)))
