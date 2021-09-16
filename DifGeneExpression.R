@@@@@@@@@@@@@@@@@@@
@ Instal packages @
@@@@@@@@@@@@@@@@@@@

@Load the CummeRbund package
@install.packages("pheatmap")
library(pheatmap)

@@@@@@@@@@@@@@@@@@@
@ Set a directory @
@@@@@@@@@@@@@@@@@@@

setwd("")

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@ Workflow for RNA sequencin data analysis @
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

@workingdir <- "Size_Factor_Estimation"
@datafiledir  <- "1_miRNA_detection_limits150k"

@datafiledir <- "R_thresholding+datastats/"
@listdir <- "1_miRNA_detection_limits150k/Thresholding"
@inputdir <- "2_normalize_5.0_95p"

myannotation="FILE1.csv"
myfile <- "mature.total.merged"
tabledir <- "Tables"
metadir <- "Metadata"
outputdir <- "DEA"


@@@@@@@@@@@@@@@@@@@@@@@@@@@
@ create output directory @
@@@@@@@@@@@@@@@@@@@@@@@@@@@

@dir.create(paste0(workingdir, sep=""),showWarnings = TRUE)
dir.create(paste0(outputdir, sep=""),showWarnings = TRUE)

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@ Loading required libraries @
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

library(sva)
@install.packages("dplyr")
library(dplyr)
@install.packages("tibble")
library(tibble)
@install.packages("data.table", repos = "https://Rdatatable.github.io/data.table", type = "source")
library(data.table)
library(DESeq2)
@install.packages("ggrepel")
library(ggrepel)
library(scales) 

dir.create(paste0(outputdir, sep=""),showWarnings = TRUE)

an <- read.csv(paste0(metadir,"/",myannotation, sep=""), header=1)
str(an)
rownames(an) <- an$ShortID
an$Batch <- as.factor(paste0("B",an$SeqBatch))
an$Adapter <- as.factor(paste0("A",an$SeqAda))
an$Index <- as.factor(an$SeqIdx)
an$Batch
an$Adapter
an$Index
SelSamples <- an 

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@ Loading merged_stats tables for decisions about processing and analysis @
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

mtables="FILE1"
mergedtables <- read.csv(paste0(metadir,"/",mtables,".csv", sep=""), row.names=1, header=1, colClasses=c("character"))
str(mergedtables)

#remove all values with commas as 1000-separators, convert all variables except "category" to numeric.
for (n in colnames(mergedtables)) {
  print(n)
  if (!(n=="category")) {
    #colnames(mergedtables) %in% n
    mergedtables[,colnames(mergedtables) == n] <- as.numeric(gsub(",","",(mergedtables[,colnames(mergedtables) == n])))
  } 
}
str(mergedtables)

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@ Creating a combined metadata and annotation dataframe @
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

mergedtables$ShortID <- as.factor(rownames(mergedtables))
mergedtables$ShortID2 <- mergedtables$ShortID

an$ShortID3 <- an$ShortID

metadata <- merge(an, mergedtables, by.x="ShortID", by.y="ShortID", all=TRUE)
metadata$ShortID2 <- factor(metadata$ShortID2, levels=unique(as.character(metadata$ShortID2)))
metadata$ShortID2 == metadata$ShortID3
metadata$ShortID2 <- NULL
metadata$ShortID3 <- NULL
mergedtables$ShortID2 <- NULL
an$ShortID3 <- NULL

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@ reating a naming scheme that R uses for colnames of data files, which includes addressing problems with leading number characters @
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

@if shortID starts with a number, let shortid4 start with an X:
@isnumber <- grepl("\\d", substr(metadata$ShortID, 1, 1))==T
@metadata$ShortID4 <- as.character(metadata$ShortID)
@metadata$ShortID4
@metadata[isnumber,]$ShortID4 <- as.character(paste0("X",metadata[isnumber,]$ShortID4,sep=""))
@metadata$ShortID4
@if shortID contains a dash ("_"), replace the dash by a dot:
@metadata$ShortID4 <- gsub("-", ".",metadata$ShortID4)
@metadata$ShortID4

SelSamples <- metadata
SelSamples$Cellline

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@ Defining, which samples to consider for processing -> for data statistics and future analyses @
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

@ remove water samples @
SelSamples <- SelSamples[!(SelSamples$CellCulture == "WATER"),]

@ consider only relevant data @
@anfilt <- data.frame(SelSamples$ShortID,SelSamples$Batch, SelSamples$Adapter, SelSamples$Index,SelSamples$Cellline,SelSamples$CellCulture,SelSamples$Replicate,SelSamples$nmiRNAperminputRNA_amolper100ng)
@anfilt <- data.frame(SelSamples$Batch, SelSamples$Adapter, SelSamples$Index,SelSamples$Techrep,SelSamples$Cellline,SelSamples$CellCulture)
anfilt <- data.frame(SelSamples$Batch, SelSamples$Adapter,SelSamples$Techrep,SelSamples$Cellline,SelSamples$CellCulture, SelSamples$miRNAp100ng.amol.,SelSamples$miRNA)

rownames(anfilt) <- SelSamples$ShortID
colnames(anfilt) <- gsub("SelSamples.","",colnames(anfilt))
colnames(anfilt) 
anfilt <- droplevels(anfilt)


temp = read.table("color_assignments/adapter.col", sep ="\t", fill=TRUE, header=FALSE, as.is=TRUE)
Adapter = temp[,2]
names(Adapter) = temp[,1]

temp = read.table("color_assignments/batch.col", sep ="\t", fill=TRUE, header=FALSE, as.is=TRUE)
Batch = temp[,2]
names(Batch) = temp[,1]

temp = read.table("color_assignments/index.col", sep ="\t", fill=TRUE, header=FALSE, as.is=TRUE)
Index = temp[,2]
names(Index) = temp[,1]

temp = read.table("color_assignments/FILE1.col", sep ="\t", fill=TRUE, header=FALSE, as.is=TRUE)
FILE1 = temp[,2]
names(FILE1) = temp[,1]
FILE1 <- FILE1[names(FILE1) %in% levels(anfilt$FILE1)]

temp = read.table("color_assignments/Cellline.col", sep ="\t", fill=TRUE, header=FALSE, as.is=TRUE)
FILE1 = temp[,2]
names(FILE1) = temp[,1]
FILE1 <- Cellline[names(FILE1) %in% levels(anfilt$FILE1)]

temp = read.table("color_assignments/Techrep.col", sep ="\t", fill=TRUE, header=FALSE, as.is=TRUE)
Techrep = temp[,2]
names(Techrep) = temp[,1]

@ now read in read counts @
@ don't use headers or row.names, since headers do not correspond to R naming conventions @
data = read.table(paste0(tabledir,"/",myfile,".txt"),stringsAsFactors = FALSE, header = FALSE, row.names = 1)

@ remove mouse entries @
data <- data[!(grepl("mmu",rownames(data))),]

@ remove rat entries @
data <- data[!(grepl("rno",rownames(data))),]

@ extract samplenames from first row of datatable @
samplenames <- as.character(data[1,])

@ extract only the string after the last occurrance of a slash @
simplenames <- sub(".*\\/", "", samplenames)
simplenames <- sub("human_and_virus_annoDB_jul2014rno-miR-325-5p", "", simplenames)
simplenames <- sub("human_and_virus_annoDB_jul2014", "", simplenames)

@ use all rows but the first, since the complicated names w. filepaths are not needed @
data <- data[2:nrow(data),]

@ apply the newly extracted simplenames @
colnames(data) <- simplenames

@ convert entries originally read in as character values into numeric values @
data1 <- as.data.frame(sapply(data, as.numeric))

@ copy over rownames which are otherwise lost @
row.names(data1) <- row.names(data) 

@ calculate sums of columns @
@colsums <- as.data.frame(t(apply(data1, 2, sum)))
@row.names(colsums) <- c("Total_counts")
@data2 <- rbind(data1, colsums) 
datafreq <- as.data.frame(apply(data1,2,function(x){x/sum(x)}))
@freqsums <- as.data.frame(t(apply(datafreq, 2, sum)))
write.csv(data1,paste0(tabledir,"/",myfile,"DEA.csv"))
write.csv(datafreq,paste0(tabledir,"/",myfile,"DEA.fr.csv"))

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@ filter out only this samples previously selected for processing @
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

datasel <- data1[,colnames(data1) %in% rownames(anfilt)]
datasel <- datasel[, order(match(colnames(datasel),rownames(anfilt)))]

freqsel <- datafreq[,colnames(datafreq) %in% rownames(anfilt)]
freqsel <- freqsel[, order(match(colnames(freqsel),rownames(anfilt)))]

@str(mat)

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@ Simplify names of clusters/merged/miRs to only numbers and essential elements @
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

@rownames(datasel) <- gsub('cluster-', '', rownames(datasel))
rownames(datasel) <- gsub('hsa-', '', rownames(datasel))
@rownames(datasel) <- gsub('mir-', '', rownames(datasel))
@rownames(datasel) <- gsub('miR-', '', rownames(datasel))
@rownames(datasel) <- gsub('STAR', '*', rownames(datasel))
@rownames(datasel) <- gsub("\\s*\\(1)","",as.character(rownames(datasel)))
remove="all others "
datasel <- datasel[!rownames(datasel) %in% remove,] 
row.names(freqsel) <- rownames(datasel)

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@ Writing out some files to document the processing and sample selection (for later data deposition) @
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

write.csv(datasel,paste0(tabledir,"/",myfile,"DEA.sel.csv",sep=""))
write.csv(freqsel,paste0(tabledir,"/",myfile,"DEA.fr.sel_selsamples.csv",sep=""))
write.csv(anfilt,paste0(metadir,"/",myfile,"DEA.meta.sel.csv",sep=""))

@@@@@@@@@@@@
@ new part @
@@@@@@@@@@@@

raw_counts = datasel

raw_frequencies = freqsel

@cali_counts = read.table(paste0(mycalfile,".txt"), sep = "\t", stringsAsFactors = FALSE, row.names = 1, 
@                         header=TRUE)

@cali_frequencies = read.table(paste0(mycalffile,".txt"), sep = "\t", stringsAsFactors = FALSE, row.names = 1, 
@                              header=TRUE)

colnames(raw_counts)
rownames(anfilt)

AnaSamples <- anfilt[row.names(anfilt) %in% colnames(raw_counts),]

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@ only use the previously selected samples @
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

@cali_counts <- cali_counts[,colnames(cali_counts) %in% AnaSamples$Samplename4]
@cali_counts <- cali_counts[, order(match(colnames(cali_counts),AnaSamples$Samplename4))]
@colnames(cali_counts) <- AnaSamples[colnames(cali_counts) == Samples4thisround$Samplename4,]$shortid

@cali_frequencies <- cali_frequencies[,colnames(cali_frequencies) %in% AnaSamples$Samplename4]
@cali_frequencies <- cali_frequencies[, order(match(colnames(cali_frequencies),AnaSamples$Samplename4)),]
@colnames(cali_frequencies) <- AnaSamples[colnames(cali_frequencies) == AnaSamples$Samplename4,]$shortid

raw_counts <- raw_counts[,colnames(raw_counts) %in% row.names(anfilt)]
raw_counts <- raw_counts[, order(match(colnames(raw_counts),row.names(anfilt))),]

raw_frequencies <- raw_frequencies[,colnames(raw_frequencies) %in% row.names(anfilt)]
raw_frequencies <- raw_frequencies[, order(match(colnames(raw_frequencies), row.names(anfilt))),]

set2 <- c("cali_07_rc","cali_11_rc","cali_12_rc","cali_14_rc","cali_15_rc","cali_16_rc","cali_26_rc","cali_28_rc","cali_31_rc","cali_35_rc")
set1 <- c("cali_01_rc","cali_04_rc","cali_17_rc","cali_18_rc","cali_20_rc","cali_24_rc","cali_25_rc","cali_27_rc","cali_43_rc","cali_44_rc")

@set2_counts <- cali_counts[set2,]
@set2_frequencies <- cali_frequencies[set2,]  
@set1_counts <- cali_counts[set1,]
@set1_frequencies <- cali_frequencies[set1,]  

row.names(anfilt) == colnames(raw_frequencies)
colnames(raw_frequencies)
rownames(raw_frequencies)

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# only consider human miRNAs used in Size factor estimation, and subsequently in count normalization @
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

@ first separate human miRNA counts from other miRNA counts @

human <- (grepl("hsa",rownames(raw_counts)))
human_counts <- raw_counts[human,]
human_frequencies <- raw_frequencies[human,]
human_counts <- raw_counts
human_frequencies <- raw_frequencies
@other_counts <- raw_counts[!human,]
@other_frequencies <- raw_frequencies[!human,]

@outputmiRNAfile=paste(workingdir,"/",outputdir,"/",tablefile,"/",tablefile,"_noc","p",percentage,experiment,sep="")
outputmiRNAfile=paste(outputdir,"/",myfile,"_noc",sep="")
outputname=paste(myfile,"_noc",sep="")

@ Idential to the keep/rowMaxs statment, but far more complicated. Also, the lower theshold is not correctly implemented. TO select 10 or more counts the value 9 has to be used @
@ install.packages("genefilter") 

@ this package is needed for the filter function pOverA; also checkout kOverA @
library(genefilter) 
f1 <- pOverA(.10, 10)

@ essentially maintain everything that has at least 10 counts in 10% of the samples. This may, however, lead to false positive results. This is just a measure to maintain detectable counts @

ffun <- filterfun(f1)
wh1 <- genefilter(raw_counts, ffun)
write.csv(wh1, file=paste0(outputmiRNAfile,"_dt10ct10p.csv"))

@ use THIS miRNA selection for subsequent normalizations @
filtered_counts <- raw_counts[wh1,]

@ Load size factors from external file @

@sizef <- read.csv(paste0(inputmiRNAfile,"_sf.csv"), row.names=1, header=1)
@sizef <- sizef [rownames(sizef) %in% SelectedSamples, ]

@v1 <- sizef
#rownames(v1) <- NULL  
@v1 <- unlist(v1)
@relist(replace(v1, v1>0, 1), skeleton=v1)
@sizef1 <- relist(replace(v1, v1>0, 1), skeleton=v1)  

@ use corrected set 1 count values instead @

combinedfilteredcounts <- rbind(filtered_counts)
@combinedfilteredcounts <- rbind(set1_counts, filtered_counts)

@ load common miRNA file @

@commonmiRNAfile=paste(workingdir,"/",listdir,"/",tablefile,"-",percentage,"0_50k.lst",sep="")
@commonmiRNAs <- read.delim(commonmiRNAfile, header=FALSE)
@commonmiRNAs <- commonmiRNAs$V1
@commonmiRNAfile="tables4R/miRNAs4scaling-1M-reads(6-6.5).txt"

@    calibrators <- grepl('cali', rownames(combinedfilteredcounts))

@  referencematrix <- rownames(filtered_counts) %in% as.vector(commonmiRNAs)
@  referencematrix <- calibrators

@tcn <- rownames(filtered_counts) %in% commonmiRNAs

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@ define DESeq parameters countData and colData @
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

@colData <- data.frame(patientdata, patient, gender, mayoclass2, samplecath, exp)

colData <- anfilt
dds_norm <- DESeqDataSetFromMatrix(countData = round(filtered_counts), 
                                   colData = colData,
                                   design = ~Cellline )

@                                design = ~gender + batchgroup)
@                                design = ~mayoclass2+gender+exp )

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@ do not use external size factors @
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

@sizeFactors(dds_norm) <- sizef1

@ do not use referencematrix for size factor estimation. Use common miRNAs detectable in all samples to achieve this @
@dds_norm <- estimateSizeFactors(dds_norm,controlGenes=tcn)

@ this will take some time, depending on the number of samples and retained miRNA @
dds_norm <- DESeq(dds_norm, parallel=TRUE, fitType="parametric", minReplicatesForReplace = 2000)

normcounts <- counts(dds_norm, normalized=TRUE)
write.csv(normcounts, file=paste0(outputmiRNAfile,"_nc.csv"))

plotDispEsts( dds_norm, ylim = c(1e-2, 1e1) )
pdf(paste0(outputmiRNAfile,"_norm.pdf"),width=5.25,height=6.0, useDingbats=FALSE)
x=rnorm(100)
y=rnorm(100,5,1)
plotDispEsts( dds_norm, ylim = c(1e-2, 1e1) )
dev.off()

@ this just shows you what are the coefficient names in the model @
resultsNames(dds_norm) 
@ the function 'results' extracts the results of interest from the dds object

@--- compar. between experiments, same groups

res.batchgroup <- results(dds_norm, contrast=c("FILE1", "CELL1", "CELL2"))  
resOrdered <- res.batchgroup[order(res.batchgroup$padj),]
head(resOrdered)
write.csv(resOrdered, file=paste0(outputmiRNAfile,"_CELL1_vs_CELL2.csv"))


@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@ Using ReportingTools to produce fancy box plots @
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

@source("http://bioconductor.org/biocLite.R")
@biocLite("ReportingTools")

library(ReportingTools)


destination <- HTMLReport(shortName = paste(outputname,"CELL1_vs_CELL2"),
                          title = paste(outputname,"miRNA that are differentially expressed betw. CELL1 vs. CELL2"),
                          reportDirectory = paste(outputmiRNAfile,"_CELL1_vs_CELL2", sep=""))

publish(res.batchgroup, destination, pvalueCutoff=1.0, DataSet=dds_norm,
        factor = colData(dds_norm)$FILE1,
        reportDir = paste(outputmiRNAfile,"_CELL1_vs_CELL2", sep=""))


@ you will find an .html file in the above path and within it you can click on the horizontal bar plots to view full vertical plots @
finish(destination) 

@@@@@@@@@@@@@@@@
@ end new part @
@@@@@@@@@@@@@@@@




