#!/usr/bin/env Rscript
##### Load require packages ##### 	
library(reshape2)
library(ggplot2)
library(VennDiagram)
############ Author Ankur Ganveer#########
##### Filename as arguments  #####
args<-commandArgs(trailingOnly=TRUE)
a<-args[1]
b<-args[2]
#####################################################
##### Check and load data for histogram of gene #####
#####################################################
size1<-file.info("overlap_gene.gtf")$size
if(size1!=0){
datag<-read.table("overlap_gene.gtf")
dfg<-data.frame(datag)}

###########################################################
##### Check and load data for histogram of transcript #####
###########################################################
size2<-file.info("overlap_trans.gtf")$size
if(size2!=0){
datat<-read.table("overlap_trans.gtf")
dft<-data.frame(datat)}

#####################################################
##### Check and load data for histogram of exon #####
#####################################################
size3<-file.info("overlap_exon.gtf")$size
if(size3!=0){
datae<-read.table("overlap_exon.gtf")
dfe<-data.frame(datae)}

################################################################
##### Load data for count of different features in <fileA> #####
################################################################
dataF1<-read.table("gtf1_count.txt", header = FALSE)
f1<-nrow(dataF1)
f1count<-dataF1[1,2] ##### Total entries in <fileA>
dataF1.df<-data.frame(dataF1[2:f1,])
names1<-list()
names1<-dataF1.df$V1
pct1<-round(dataF1.df$V2/sum(dataF1.df$V2)*100,digits=2)
names1<-paste(names1,"\n",pct1)
names1<-paste(names1,"%",sep="")

################################################################
##### Load data for count of different features in <fileB> #####
################################################################
dataF2<-read.table("gtf2_count.txt", header = FALSE)
f2<-nrow(dataF2)
f2count<-dataF2[1,2] ##### Total entries in <fileB>
dataF2.df<-data.frame(dataF2[2:f2,])
names2<-dataF2.df$V1
pct2<-round(dataF2.df$V2/sum(dataF2.df$V2)*100,digits=2)
names2<-paste(names2,"\n",pct2)
names2<-paste(names2,"%",sep="")

################################################################
##### Load data for count of overlapped different features #####
################################################################
dataOvr<-read.table("overlap_count.txt",header = FALSE)
Ovr<-nrow(dataOvr)
ovrcount<-dataOvr[1,2] ##### Total overlapped features
dataOvr.df<-data.frame(dataOvr[2:Ovr,])
namesOvr<-dataOvr.df$V1
pctOvr<-round(dataOvr.df$V2/sum(dataOvr.df$V2)*100,digits=2)
namesOvr<-paste(namesOvr,"\n",pctOvr)
namesOvr<-paste(namesOvr,"%",sep="")

####################################################################
##### Load data for count of non-overlapped different features #####
####################################################################
dataNovr<-read.table("noverlap_count.txt", header = FALSE)
Novr<-nrow(dataNovr)
novrcount<-dataNovr[1,2] ##### Total non-overlapped features
dataNovr.df<-data.frame(dataNovr[2:Novr,])
namesNovr<-dataNovr.df$V1
pctNovr<-round(dataNovr.df$V2/sum(dataNovr.df$V2)*100,digits=2)
namesNovr<-paste(namesNovr,"\n",pctNovr)
namesNovr<-paste(namesNovr,"%",sep="")

####################################################################
##### Mergering the counts of different features from <fileA>  #####
##### and <fileB> with the count of overlapped different       #####
##### features in one data set for creating VennDiagrams.      #####
####################################################################
dataV.df<-data.frame(dataF1)
dataV.df$V3<-dataF2.df$V2[match(dataV.df$V1,dataF2.df$V1)]
dataV.df$V4<-dataOvr.df$V2[match(dataV.df$V1,dataOvr.df$V1)]
dataV1.df<-dataV.df[complete.cases(dataV.df),]
n<-nrow(dataV1.df)

####################################################################
##### Check and swap the count of a feature in crossover area, ##### 
##### when it is greater than minimum of the count of same     #####
##### feaure from 2 input files, with the minimum count of     #####
##### same feature from either file. 			       #####
####################################################################
for (i in 1:n){
if(dataV1.df$V2[i]<dataV1.df$V3[i]){
min<-dataV1.df$V2[i]
}
else
min<-dataV1.df$V3[i]
if((dataV1.df$V4[i] > dataV1.df$V2[i]) || (dataV1.df$V4[i] > dataV1.df$V3[i])){
dataV1.df$V4[i]<-min}
}

####################################################################
##### Creating dataset for total entries barplot 	       #####
####################################################################
ovrcount1<-sum(dataV1.df$V4)
filename<-c(a,b,'overlap','noverlap')
linecount<-c(f1count,f2count,ovrcount1,novrcount)
filecount.df<-data.frame(filename,linecount)

dataF1B.df<-dataF1.df
dataF2B.df<-dataF2.df
dataOvrB.df<-dataOvr.df
dataNovrB.df<-dataNovr.df

####################################################################
##### Generating names of 4 bars for total entries bar plot	   #####
####################################################################
colnames(dataF1B.df)[2]<-a ##### name of <fileA>
colnames(dataF2B.df)[2]<-b ##### name of <fileB>
colnames(dataOvrB.df)[2]<-"Overlap"
colnames(dataNovrB.df)[2]<-"No-Overlap"

####################################################################
##### Creating dataset for counts of different features        ##### 
##### of <fileA> <fileB> barplot                               #####
####################################################################
mer<-merge(dataF1B.df,dataF2B.df,by="V1")
df3<-melt(mer)
colnames(df3)<-c("Features", "Files", "Count")

####################################################################
##### Creating dataset for counts of different features        ##### 
##### for overlapped and non-overlapped  barplot               #####
####################################################################
mer1<-merge(dataOvrB.df,dataNovrB.df,by="V1")
df4<-melt(mer1)
colnames(df4)<-c("Features", "Condition", "Count")

####################################################################
##### Creating dataset for counts of different features        ##### 
##### breakdown by chromosomes of <fileA> for barplot          #####
####################################################################
dataChrF1<-read.csv("chr_gtf1.txt", header=FALSE)
tmp1<-split(dataChrF1,dataChrF1$V2)
dataChrF1N<-lapply(1:length(tmp1),function(x) assign(paste("dataChrF1.",x,sep=""),tmp1[[x]],envir=.GlobalEnv))

####################################################################
##### Creating dataset for counts of different features        ##### 
##### breakdown by chromosomes of <fileB> for barplot          #####
####################################################################
dataChrF2<-read.csv("chr_gtf2.txt", header=FALSE)
tmp2<-split(dataChrF2,dataChrF2$V2)
dataChrF2N<-lapply(1:length(tmp2),function(x) assign(paste("dataChrF2.",x,sep=""),tmp2[[x]],envir=.GlobalEnv))

####################################################################
##### Creating dataset for counts of different features        ##### 
##### breakdown by strands of <fileA> for barplot              #####
####################################################################
dataSndF1<-read.csv("snd_gtf1.txt", header=FALSE)
tmp3<-split(dataSndF1,dataSndF1$V2)
dataSndF1N<-lapply(1:length(tmp3),function(x) assign(paste("dataSndF1.",x,sep=""),tmp3[[x]],envir=.GlobalEnv))

####################################################################
##### Creating dataset for counts of different features        ##### 
##### breakdown by strands of <fileB> for barplot              #####
####################################################################
dataSndF2<-read.csv("snd_gtf2.txt", header=FALSE)
tmp4<-split(dataSndF2,dataSndF2$V2)
dataSndF2N<-lapply(1:length(tmp4),function(x) assign(paste("dataSndF2.",x,sep=""),tmp4[[x]],envir=.GlobalEnv))

####################################################################
##### Creating dataset for counts of different overlapped      ##### 
##### features breakdown by chromosomes for barplot            #####
####################################################################
dataChrOvr<-read.csv("chr_file2.txt", header=FALSE)
tmp5<-split(dataChrOvr,dataChrOvr$V2)
dataChrOvrN<-lapply(1:length(tmp5),function(x) assign(paste("dataChrOvr.",x,sep=""),tmp5[[x]],envir=.GlobalEnv))

####################################################################
##### Creating dataset for counts of different non overlapped  ##### 
##### features breakdown by chromosomes for barplot            #####
####################################################################
dataChrNovr<-read.csv("chr_file1.txt", header=FALSE)
tmp6<-split(dataChrNovr,dataChrNovr$V2)
dataChrNovrN<-lapply(1:length(tmp6),function(x) assign(paste("dataChrNovr.",x,sep=""),tmp6[[x]],envir=.GlobalEnv))

####################################################################
##### Creating dataset for counts of different overlapped      ##### 
##### features breakdown by strands for barplot                #####
####################################################################
dataSndOvr<-read.csv("snd_file2.txt", header=FALSE)
tmp7<-split(dataSndOvr,dataSndOvr$V2)
dataSndOvrN<-lapply(1:length(tmp7),function(x) assign(paste("dataSndOvr.",x,sep=""),tmp7[[x]],envir=.GlobalEnv))

####################################################################
##### Creating dataset for counts of different non overlapped  ##### 
##### features breakdown by strands for barplot                #####
####################################################################
dataSndNovr<-read.csv("snd_file1.txt", header=FALSE)
tmp8<-split(dataSndNovr,dataSndNovr$V2)
dataSndNovrN<-lapply(1:length(tmp8),function(x) assign(paste("dataSndNovr.",x,sep=""),tmp8[[x]],envir=.GlobalEnv))

##### Title for piechart #####
pietitle1<-paste(a,"- Feature Distribution")
pietitle2<-paste(b,"- Feature Distribution")


####################################################################
##### Creating dataset for chromosome venndiagram              ##### 
####################################################################
chr1<-read.table("chr_intern1.txt")
chr2<-read.table("chr_intern2.txt")
chr1<-chr1[,c(2,1)]
chr2<-chr2[,c(2,1)]
cchr1<-nrow(chr1)
cchr2<-nrow(chr2)
data<-chr1
data$V3<-chr2$V1[match(data$V2,chr2$V2)] 
data1<-data[complete.cases(data),]  
cchr<-nrow(data1)
colnames(chr1)[2]<-a
colnames(chr2)[2]<-b

merchr<-merge(chr1,chr2,by="V2")
chrmelt<-melt(merchr)
colnames(chrmelt)<-c("Chromosomes", "Files", "Count")

####################################################################
#####              Creating first PDF file		       ##### 
####################################################################
pdf("Plots.pdf")

##### Bar plot for total entries in <fileA>, <fileB>, and found as overlap and non-overlap #####
ggplot(filecount.df,aes(x=filename,y=linecount))+geom_bar(stat="identity",position="dodge")+ ggtitle("Total Entries")+xlab("Input files and condition")+ylab("Number of Entries")

##### Bar plot for total entries by features in <fileA> and <fileB> #####
ggplot(df3,aes(Features,Count,fill=Files))+geom_bar(stat="identity",position="dodge")+ggtitle("Counts of Different Features")
##### Bar plot for total entries by features found as overlapped and non-overlapped #####
ggplot(df4,aes(Features,Count,fill=Condition))+geom_bar(stat="identity",position="dodge")+ggtitle("Counts of Different Features")

##### Pie diagrams #####
pie(dataF1.df$V2,labels=names1,col=rainbow(length(names1)),main=pietitle1)
pie(dataF2.df$V2,labels=names2,col=rainbow(length(names2)),main=pietitle2)
#pie(dataOvr.df$V2,labels=namesOvr,col=rainbow(length(namesOvr)),main="Feature Distribution of Total Overlaps Found")
#pie(dataNovr.df$V2,labels=namesNovr,col=rainbow(length(namesNovr)),main="Feature Distribution of Total Non-Overlaps Found")

##### Venn Diagrams for each different features #####
for(i in 1:n){
grid.newpage()
write(dataV1.df$V1[i])
venn.plot <- draw.pairwise.venn(
area1 = dataV1.df$V2[i],
area2 = dataV1.df$V3[i],
cross.area = dataV1.df$V4[i],
category = c(a, b),
fill = c("blue", "red"),
lty = "blank",
cex = 2,
cat.cex = 2,
cat.pos = c(285, 105),
cat.dist = 0.09,
cat.just = list(c(-1, -1), c(1, 1)),
ext.pos = 30,
ext.dist = -0.05,
ext.length = 0.85,
ext.line.lwd = 2,
ext.line.lty = "dashed"
);
grid.draw(venn.plot)
grid.text(dataV1.df$V1[i])
}

##### Venn Diagram for chromosome #####
grid.newpage()
venn.plot<-draw.pairwise.venn(cchr1,cchr2,cchr,category = c(a, b),fill = c("blue", "red"),lty = "blank");
grid.text("Chromosomes")
grid.draw(venn.plot)

##### Bar plot for counts of different common chromosomes in <fileA> <fileB> #####
ggplot(chrmelt,aes(Chromosomes,Count,fill=Files))+geom_bar(stat="identity",position="dodge")+ ggtitle("Counts of Common Chromosomes")+theme(axis.text.x = element_text(angle = 90, hjust = 1))

##### Ggplot Histogram - GENES #####
if(size1!=0){
print(qplot(dfg$V14, geom="histogram",  main="Histogram of Differences in start positions of overlapping GENES", xlab="Difference in Start Positions"))
print(qplot(dfg$V15, geom="histogram", main="Histogram of Differences in end positions of overlapping GENES", xlab="Difference in End Positions"))
}

##### Ggplot Histogram - TRANSCRIPTS#####
if(size2!=0){
print(qplot(dft$V14, geom="histogram",main="Histogram of Differences in Start Position of overlapping TRANSCRIPTS", xlab="Difference in Start Positions"))
print(qplot(dft$V15, geom="histogram",main="Histogram of Differences in End Position of overlapping TRANSCRIPTS", xlab="Difference in End Positions"))
}

##### Ggplot Histogram - EXONS #####
if(size3!=0){
print(qplot(dfe$V14, geom="histogram",binwidth=40,main="Histogram of Differences in Start Positions of overlapping EXONS", xlab="Difference in Start Positions"))
print(qplot(dfe$V15, geom="histogram",binwidth=40,main="Histogram of Differences in End Positions of overlapping EXONS", xlab="Difference in End Positions"))
}

##### Hist Histogram - GENES #####
if(size1!=0){
print(hist(dfg$V14,breaks=40,xlab="Difference in Start Positions",main="Histogram of Differences in Start Positions of overlapping GENES"));
print(hist(dfg$V15,breaks=40,xlab="Difference in End Positions",main="Histogram of Differences in End Positions of overlapping GENES"));
}
##### Hist Histogram - TRANSCRIPTS #####
if(size2!=0){
print(hist(dft$V14,breaks=40,xlab="Difference in Start Positions",main="Histogram of Differences in Start Positions of overlapping TRANSCRIPTS"));
print(hist(dft$V15,breaks=40,xlab="Difference in End Positions",main="Histogram of Differences in End Positions of overlapping TRANSCRIPTS"));
}
##### Hist Histogram - EXONS #####
if(size3!=0){
print(hist(dfe$V14,breaks=40,xlab="Difference in Start Positions",main="Histogram of Differences in Start Positions of overlapping EXONS"));
print(hist(dfe$V15,breaks=40,xlab="Difference in End Positions",main="Histogram of Differences in End Positions of overlapping EXONS"));
}

#### Bar plot for counts of features breakdown by chromosome #####
print(ggplot(dataChrF1,aes(V2,V3,fill=V1))+ geom_bar(stat="identity",position="dodge")+xlab ("Different Features") +ylab("Count of Features") + ggtitle(paste("Counts of Features breakdown by Chromosomes in the file-",a))+theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_discrete(name = "Chromosomes")) 
print(ggplot(dataChrF2,aes(V2,V3,fill=V1))+ geom_bar(stat="identity",position="dodge")+xlab ("Different Features") +ylab("Count of Features") + ggtitle(paste("Counts of Features breakdown by Chromosomes in the file-",b))+theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_discrete(name = "Chromosomes"))
#### Bar plot for counts of features breakdown by strand #####
print(ggplot(dataSndF1,aes(V2,V3,fill=V1))+ geom_bar(stat="identity",position="dodge")+xlab ("Different Features") +ylab("Count of Features") + ggtitle(paste("Counts of Features breakdown by Strands in the file-",a))+theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_discrete(name = "Strands"))
print(ggplot(dataSndF2,aes(V2,V3,fill=V1))+ geom_bar(stat="identity",position="dodge")+xlab ("Different Features") +ylab("Count of Features") + ggtitle(paste("Counts of Features breakdown by Strands in the file-",b))+theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_discrete(name = "Strands"))
#### Bar plot for counts of features breakdown by chromosome #####
print(ggplot(dataChrOvr,aes(V2,V3,fill=V1))+ geom_bar(stat="identity",position="dodge")+xlab ("Different Features") +ylab("Count of Features") + ggtitle("Counts of Overlapping Features breakdown by Chromosomes")+theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_discrete(name = "Chromosomes"))
print(ggplot(dataChrNovr,aes(V2,V3,fill=V1))+ geom_bar(stat="identity",position="dodge")+xlab ("Different Features") +ylab("Count of Features") + ggtitle("Counts of Non Overlapping Features breakdown by Chromosomes")+theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_discrete(name = "Chromosomes"))
#### Bar plot for counts of features breakdown by strand #####
print(ggplot(dataSndOvr,aes(V2,V3,fill=V1))+ geom_bar(stat="identity",position="dodge")+xlab ("Different Features") +ylab("Count of Features") + ggtitle("Counts of Overlapping Features breakdown by Strands")+theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_discrete(name = "Strands"))
print(ggplot(dataSndNovr,aes(V2,V3,fill=V1))+ geom_bar(stat="identity",position="dodge")+xlab ("Different Features") +ylab("Count of Features") + ggtitle("Counts of Non Overlapping Features breakdown by Strands")+theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_discrete(name = "Strands"))

dev.off()

####################################################################
#####              Creating second DetailedPDF file	       ##### 
#################################################################### 
pdf("DetailPlots.pdf")

##### Count of features by different chromosomes for <fileA> #####
print(ggplot(dataChrF1,aes(V1,V3,fill=V2))+ geom_bar(stat="identity",position="dodge")+xlab ("Different Chromosomes") +ylab("Count of Features") + ggtitle(paste("Counts of Features breakdown by Chromosomes in the file-",a))+theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_discrete(name = "Features"))

##### Separate bar plots for Count of separate features by different chromosomes for <fileA> #####
for(i in 1:length(tmp1)){
#print(ggplot(dataChrF1N[[i]],aes(V2,V3,fill=V1))+ geom_bar(stat="identity",position="dodge")+xlab ("Feature") +ylab("Count of Features") + ggtitle(paste("Counts of Features breakdown by Chromosomes in the file-",a))+theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_discrete(name = "Chromosomes"))
print(ggplot(dataChrF1N[[i]],aes(V1,V3,fill=V2))+ geom_bar(stat="identity",position="dodge")+xlab ("Different Chromosomes") +ylab("Count of Features") + ggtitle(paste("Counts of Features by Chromosomes in the file-",a))+theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_discrete(name = "Feature"))
}

##### Count of features by different strands for <fileA> #####
print(ggplot(dataSndF1,aes(V1,V3,fill=V2))+ geom_bar(stat="identity",position="dodge")+xlab ("Different Strands") +ylab("Count of Features") + ggtitle(paste("Counts of Features breakdown by Strands in the file-",a))+theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_discrete(name = "Features"))

#for(i in 1:length(tmp3)){
#print(ggplot(dataSndF1N[[i]],aes(V2,V3,fill=V1))+ geom_bar(stat="identity",position="dodge")+xlab ("Different Features") +ylab("count of features") + ggtitle(paste("Counts of",a,"Features by Strands"))+theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_discrete(name = "Strands"))
#print(ggplot(dataSndF1N[[i]],aes(V1,V3,fill=V2))+ geom_bar(stat="identity",position="dodge")+xlab ("Different Chromosomes") +ylab("count of features") + ggtitle(paste("Counts of",a,"Strand by Features"))+theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_discrete(name = "Features"))
#}

##### Count of features by different chromosomes for <fileB> #####
print(ggplot(dataChrF2,aes(V1,V3,fill=V2))+ geom_bar(stat="identity",position="dodge")+xlab ("Different Chromosomes") +ylab("Count of Features") + ggtitle(paste("Counts of Features breakdown by Chromosomes in the file-",b))+theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_discrete(name = "Features"))

##### Separate bar plots for Count of separate features by different chromosomes for <fileB> #####
for(i in 1:length(tmp2)){
#print(ggplot(dataChrF2N[[i]],aes(V2,V3,fill=V1))+ geom_bar(stat="identity",position="dodge")+xlab ("Feature") +ylab("Count of Features") + ggtitle(paste("Counts of Feature breakdown by Chromosomes in the file-",b))+theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_discrete(name = "Chromosomes"))
print(ggplot(dataChrF2N[[i]],aes(V1,V3,fill=V2))+ geom_bar(stat="identity",position="dodge")+xlab ("Different Chromosomes") +ylab("count of features") + ggtitle(paste("Counts of Features by Chromosomes in the file-",b))+theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_discrete(name = "Feature"))
}

##### Count of features by different strands for <fileB> #####
print(ggplot(dataSndF2,aes(V1,V3,fill=V2))+ geom_bar(stat="identity",position="dodge")+xlab ("Different Strands") +ylab("Count of Features") + ggtitle(paste("Counts of Features breakdown by Strands in the file-",b))+theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_discrete(name = "Features"))

#for(i in 1:length(tmp4)){
#print(ggplot(dataSndF2N[[i]],aes(V2,V3,fill=V1))+ geom_bar(stat="identity",position="dodge")+xlab ("Different Features") +ylab("count of features") + ggtitle(paste("Counts of",b,"Features by Strands"))+theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_discrete(name = "Strands"))
#print(ggplot(dataSndF2N[[i]],aes(V1,V3,fill=V2))+ geom_bar(stat="identity",position="dodge")+xlab ("Different Chromosomes") +ylab("count of features") + ggtitle(paste("Counts of",b,"Strand by Feature"))+theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_discrete(name = "Features"))
#}

##### Count of features by different chromosomes found as overlapped #####
print(ggplot(dataChrOvr,aes(V1,V3,fill=V2))+ geom_bar(stat="identity",position="dodge")+xlab ("Different Chromosomes") +ylab("Count of Features") + ggtitle("Counts of Overlapping Features breakdown by Chromosomes")+theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_discrete(name = "Features"))

##### Separate bar plots for Count of separate features by different chromosomes found as overlapped #####
for(i in 1:length(tmp5)){
#print(ggplot(dataChrOvrN[[i]],aes(V2,V3,fill=V1))+ geom_bar(stat="identity",position="dodge")+xlab ("Feature") +ylab("Count of Features") + ggtitle("Counts of Overlapping Features breakdwon by Chromosomes")+theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_discrete(name = "Chromosomes"))
print(ggplot(dataChrOvrN[[i]],aes(V1,V3,fill=V2))+ geom_bar(stat="identity",position="dodge")+xlab ("Different Chromosomes") +ylab("count of features") + ggtitle("Counts of Overlapping Features breakdown by Chromosomes")+theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_discrete(name = "Features"))
}

##### Count of features by different strands found as overlapped #####
print(ggplot(dataSndOvr,aes(V1,V3,fill=V2))+ geom_bar(stat="identity",position="dodge")+xlab ("Different Strands") +ylab("Count of Features") + ggtitle("Counts of Overlapping Features breakdown by Strands")+theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_discrete(name = "Features"))

#for(i in 1:length(tmp7)){
#print(ggplot(dataSndOvrN[[i]],aes(V2,V3,fill=V1))+ geom_bar(stat="identity",position="dodge")+xlab ("Different Features") +ylab("count of features") + ggtitle("Counts of Overlapping Features by Strands")+theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_discrete(name = "Strands"))
#print(ggplot(dataSndOvrN[[i]],aes(V1,V3,fill=V2))+ geom_bar(stat="identity",position="dodge")+xlab ("Different Chromosomes") +ylab("count of features") + ggtitle("Counts of Overlapping Strand by Features")+theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_discrete(name = "Features"))
#}

##### Count of features by different chromosomes found as overlapped #####
print(ggplot(dataChrNovr,aes(V1,V3,fill=V2))+ geom_bar(stat="identity",position="dodge")+xlab ("Different Chromosomes") +ylab("Count of Features") + ggtitle("Counts of Non Overlapping Features breakdown by Chromosomes")+theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_discrete(name = "Features"))

##### Separate bar plots for Count of separate features by different chromosomes found as overlapped #####
for(i in 1:length(tmp6)){
#print(ggplot(dataChrNovrN[[i]],aes(V2,V3,fill=V1))+ geom_bar(stat="identity",position="dodge")+xlab ("Feature") +ylab("Count of Features") + ggtitle("Counts of Non Overlapping Features breakdown by Chromosome")+theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_discrete(name = "Chromosomes"))
print(ggplot(dataChrNovrN[[i]],aes(V1,V3,fill=V2))+ geom_bar(stat="identity",position="dodge")+xlab ("Different Chromosomes") +ylab("count of features") + ggtitle("Counts of Non Overlapping Features breakdown by Chromosomes")+theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_discrete(name = "Features"))
}

##### Count of features by different strands found as overlapped #####
print(ggplot(dataSndNovr,aes(V1,V3,fill=V2))+ geom_bar(stat="identity",position="dodge")+xlab ("Different Strands") +ylab("Count of Features") + ggtitle("Counts of Non Overlapping Features breakdown by Strands")+theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_discrete(name = "Features"))

#for(i in 1:length(tmp8)){
#print(ggplot(dataSndNovrN[[i]],aes(V2,V3,fill=V1))+ geom_bar(stat="identity",position="dodge")+xlab ("Different Features") +ylab("count of features") + ggtitle("Counts of Non Overlapping Features by Strands")+theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_discrete(name = "Strands"))
#print(ggplot(dataSndNovrN[[i]],aes(V1,V3,fill=V2))+ geom_bar(stat="identity",position="dodge")+xlab ("Different Chromosomes") +ylab("count of features") + ggtitle("Counts of Non Overlapping Strand by Features")+theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_discrete(name = "Features"))
#}

dev.off()