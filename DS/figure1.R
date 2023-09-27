library(pheatmap)
library(tidyverse)
library(readxl)
library(rlist)
library(ggplot2)
require("ggplot2")
library(ggpubr)

setwd('/research_jude/rgs01_jude/groups/xmagrp/projects/SingleCellDNA/common/DS-AML/WGS/tracking/Zhikai/figures/figure1_proportion//')
## refseq parsing ##
refdir = "/research_jude/rgs01_jude/groups/xmagrp/projects/SensitiveDetection/xmagrp/code/mRNA/ZL/reference/"
refseq38 = read.table(paste(refdir,"ncbiRefSeqCurated_hg38_default.txt",sep=""),sep="\t",head=TRUE)
mygene = refseq38 %>% filter(gName=="GATA1",isDefault=="Yes") %>% select(exonS, exonE)
exonSlst = mygene$exonS
exonElst = mygene$exonE
##
buffer = 5
E2start = as.integer(strsplit(exonSlst,",")[[1]][2])+1-buffer
E2end = as.integer(strsplit(exonElst,",")[[1]][2])+buffer
##
diff = E2end-E2start+1
allpos = data.frame(c(E2start:E2end))
colnames(allpos) = c("stpos")
##
## markers parsing ##
input = read_excel("../DS_AML_zl.xlsx",sheet = "markers")
gata1 = input %>% select(Marker_hg38, Type, Gene) %>% filter(Gene == "GATA1")
ITD = gata1 %>% filter(grepl("ITD",Type)) %>% 
  separate(Marker_hg38, into=c("chr","pos","ref","alt"), extra = 'drop', remove = FALSE) %>% 
  mutate(stpos = as.integer(pos)-1-nchar(alt)+2) %>% 
  mutate(sppos = as.integer(pos)-1)  %>% rownames_to_column(var="yy") %>%
  mutate(newType = "ITD") %>% select(newType, stpos, sppos, yy)
IndelDel = gata1 %>% filter(Type == "IndelDel") %>% 
   separate(Marker_hg38, into=c("chr","pos","ref","alt"), extra = 'drop', remove = FALSE) %>% 
  mutate(stpos = as.integer(pos),
         sppos = as.integer(pos)+nchar(ref)) %>% 
  rownames_to_column(var="yy") %>%
  mutate(newType = Type) %>% select(newType, stpos, sppos, yy)
IndelIns = gata1 %>% filter(Type == "IndelIns") %>% 
  separate(Marker_hg38, into=c("chr","pos","ref","alt"), extra = 'drop', remove = FALSE) %>% 
  mutate(stpos = as.integer(pos), 
         sppos = as.integer(pos)+nchar(ref)) %>% rownames_to_column(var="yy") %>%
  mutate(newType = Type) %>% select(newType, stpos, sppos, yy)
IndelComplex = gata1 %>% filter(Type == "IndelComplex") %>% 
  separate(Marker_hg38, into=c("chr","pos","ref","alt"), extra = 'drop', remove = FALSE) %>% 
  mutate(stpos = as.integer(pos), sppos = as.integer(pos)+nchar(ref)) %>% 
  rownames_to_column(var="yy") %>%
  mutate(newType = Type) %>% select(newType, stpos, sppos, yy)
SNV = gata1 %>% filter(Type == "SNV" | Type == "MNV") %>% 
  separate(Marker_hg38, into=c("chr","pos","ref","alt"), extra = 'drop', remove = FALSE) %>% 
  mutate(stpos = as.integer(pos), sppos = as.integer(pos)+nchar(ref)) %>% 
  rownames_to_column(var="yy") %>%
  mutate(newType = ifelse(Type == "MNV", "SNV", Type)) %>% select(newType, stpos, sppos, yy)
SV = gata1 %>% filter(Type == "SV") %>% 
  separate(Marker_hg38, sep="\\." , into=c("chrA","posA","strA","chrB","posB","strB"), 
           extra = 'drop', remove = FALSE) %>% 
  mutate(stpos = as.integer(posA)) %>% mutate(sppos = as.integer(posB)) %>% rownames_to_column(var="yy") %>%
  mutate(newType = "SV") %>% select(newType, stpos, sppos, yy) 
#alltypes = rbind(ITD, IndelDel, IndelIns, IndelComplex, SNV, SV)
alltypes = rbind(ITD, IndelDel, IndelIns, IndelComplex, SNV)
alltypes$yy = as.integer(as.character(alltypes$yy))

scaleFUN = function(x) sprintf("%.1f", x)
p1 = ggplot(data=alltypes, aes(color=newType))+
  geom_segment(aes(x = stpos, y = yy, xend = sppos, yend = yy), 
               linewidth = 1)+scale_y_continuous(labels=scaleFUN)+
  coord_cartesian(xlim=c(E2start, E2end))+
  theme_classic(base_size = 18)+xlab("")+
  theme(axis.text.x=element_blank(),
       axis.ticks.x=element_blank(),
       axis.title.x=element_blank(),
       plot.margin=unit(c(0,0,0,0),"cm"))+
  geom_vline(xintercept = c(E2start+buffer, E2end-buffer), 
            color = "blue", linewidth=1)+ylab("Events")
p1
#####
totalevent = nrow(alltypes1)
alltypes1 = alltypes %>% filter(newType != "SV") %>% select(newType, stpos)
alltypes2 = alltypes1 %>% group_by(stpos, newType) %>% summarise(n=n()) %>% ungroup()
#
diff = E2end-E2start+1
allpos = data.frame(c(E2start:E2end))
colnames(allpos) = c("stpos")
binsize = 5
allpos = data.frame(seq(E2start, E2end+binsize, binsize))
colnames(allpos) = "stpos"
#
typebase = c("IndelDel","SNV","IndelIns","ITD","IndelComplex")
newst = c()
newsp = c()
newtype = c()
newnum = c()
for (r in c(1:(nrow(allpos)-1))){
  st = allpos[r, "stpos"]
  sp = allpos[r+1, "stpos"]
  print (paste(st, sp, sep = ","))
  del = 0
  ins = 0
  itd = 0
  snv = 0
  com = 0
  for (rr in c(1:nrow(alltype2))){
    type = as.character(alltype2[rr, "newType"])
    number = alltype2[rr, "n"]
    pos = alltype2[rr, "stpos"]
    if (pos >= st & pos < sp){
      if (type == "IndelDel"){
        del = del + number
      } else if (type == "SNV"){
        snv = snv + number
      } else if (type == "IndelIns"){
        ins = ins + number
      } else if (type == "ITD"){
        itd = itd + number
      } else if (type == "IndelComplex"){
        com = com + number
      }
    }
  }
  for (tt in typebase){
    newst = c(newst, st)
    newsp = c(newsp, sp)
    if (tt == "IndelDel"){
      newtype = c(newtype, tt)
      newnum = c(newnum, del)
    } else if (tt == "SNV"){
      newtype = c(newtype, tt)
      newnum = c(newnum, snv)
    } else if (tt == "IndelIns"){
      newtype = c(newtype, tt)
      newnum = c(newnum, ins)
    } else if (tt == "ITD"){
      newtype = c(newtype, tt)
      newnum = c(newnum, itd)
    } else if (tt == "IndelComplex"){
      newtype = c(newtype, tt)
      newnum = c(newnum, com)
    }
  }
}
final = data.frame(cbind(newst, newsp, newtype, newnum))
head(final)
final$newst = as.integer(as.character(final$newst))
final$newsp = as.integer(as.character(final$newsp))
final$newnum = as.integer(as.character(final$newnum))
rownames(final) = seq(1:nrow(final))
final1 = final %>% arrange(by=newst)
ggplot(final, aes(x=newst, y=newnum, group=newtype)) +
  geom_line(aes(color=newtype))

alltypes3 = allpos %>% left_join(alltypes2, by="stpos")
final %>% arrange(newst, newtype)
#######
fg = alltypes1 %>% group_by(newType, stpos) %>% summarise(total=n()) %>% ungroup()
####
fg1 = data.frame()
for (i in unique(fg$newType)){
  x = fg %>% filter(newType == i)
  r = allpos %>% left_join(x, by="stpos") %>% 
    mutate(value = ifelse(is.na(total),0,total),
           type = i) %>% select(stpos, value, type)
  fg1 = rbind(fg1, r)
}
fg2 = fg1 %>% group_by(type) %>% mutate(cumtype = cumsum(value)) %>% ungroup()
bg = alltypes1 %>% group_by(stpos) %>% summarise(total=n())
bg1 = allpos %>% left_join(bg, by="stpos") %>% 
  mutate(value = ifelse(is.na(total),0,total)) %>% select(-total) %>% 
  mutate(cumy=cumsum(value))
df1 = fg2 %>% left_join(bg1, by="stpos") %>% mutate(tmp = cumtype/cumy) %>%
  mutate(prop = ifelse(is.na(tmp),0,tmp)) %>% select(-tmp)

p2 = ggplot(data=df1,aes(x=stpos, y=prop, group=type)) +
  geom_line(aes(color=type))+
  theme_classic(base_size = 18) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
  xlim(E2start,E2end)+theme(plot.margin=unit(c(0,0,0,0),"cm"))+
  geom_vline(xintercept = c(E2start+buffer, E2end-buffer), 
             color = "blue", size=1)+ylab("Event enrichment")+
  xlab("Position")
ggarrange(p1, p2, nrow=2)

