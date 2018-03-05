library(openxlsx)
library(plyr)
library(dplyr)
library(data.table)
library(caret)
library(cvTools)

## get the job ID
r <- as.numeric(Sys.getenv("SGE_TASK_ID"))	
# r is the interger from 1 to 10 

# the number of CV folds
cvRuns = 10

# load SNPs
load("SNPs_20171024.rda")
#==============
# input files
# file that contains information about blacklist mutations to exclude
#blacklistFile = "blacklist_from_July-17-2017.xlsx"
blacklistFile = "blacklist_from_October-23-2017.xlsx"

# file with WBC data
wbcFile = "WBC Mutations (Filtered + Optical Duplicates Removed + Q30), October 6, 2017.csv"
# file with plasmaSeq data to analyze
#plasmaSeqDataFile = "MEGA PlasmaSeq Data for Lu, Filtered Mutations (Optical Duplicates Removed + Q30) with PT and WBC Data, October 6, 2017.csv"
plasmaSeqDataFile = "MEGA PlasmaSeq Data for Lu, Filtered Mutations (Optical Duplicates Removed + Q30) with PT and WBC Data, December 8, 2017.csv"

#=================================
# START ANALYSIS
# blacklist mutations
blacklist = read.xlsx(blacklistFile, sheet = 1)
blacklist$mutID = paste(blacklist$Chrom,blacklist$Position,blacklist$BaseFrom,blacklist$BaseTo)
blacklist = blacklist$mutID

#==============
# WBC data from separate file
#==============
wbcDataFromFile = fread(wbcFile)
wbcDataFromFile = data.frame(wbcDataFromFile)
wbcDataFromFile = wbcDataFromFile[!(wbcDataFromFile$ampMatchName %in% amplicon),]
colnames(wbcDataFromFile)[c(20,14,17)] = c("TotalUID", "UID1","UID2")
## all plasma samples with matched WBC
matched = unique(wbcDataFromFile$Matched.Plasma)
matched = matched[!is.na(matched)]

wbcDataFromFile$mutID = paste(wbcDataFromFile$Chrom,wbcDataFromFile$Position,wbcDataFromFile$BaseFrom,wbcDataFromFile$BaseTo)

## deleting NA rows and snps from WBC mutations
wbcDataFromFile = wbcDataFromFile[!is.na(wbcDataFromFile$Position),]
wbcDataFromFile = wbcDataFromFile[!(wbcDataFromFile$mutID %in% snps),]

## remove duplicated row (same mutation in the same template)
wbcDataFromFile$remove = paste(wbcDataFromFile$Template,wbcDataFromFile$mutID)

dupPairs <- names(which(table(wbcDataFromFile$remove)>1))
newuniquenewnormal <- wbcDataFromFile%>% 
  filter(remove %in% dupPairs)%>%
  group_by(remove)%>%
  slice(which.max(TotalUID)) %>% 
  as.data.frame()

wbcDataFromFile <- rbind(newuniquenewnormal,
                         wbcDataFromFile%>% 
                           filter(!(remove %in% dupPairs)))

wbcDataFromFile = data.frame(wbcDataFromFile)

## numeric columns
col = c("SM1","UID1","SM2","UID2","TotalSM","TotalUID")

wbcDataFromFile[,col] = apply(wbcDataFromFile[,col], 2, function(x) as.numeric(x))
wbcDataFromFile$maf1 = wbcDataFromFile$SM1/wbcDataFromFile$UID1
wbcDataFromFile$maf2 = wbcDataFromFile$SM2/wbcDataFromFile$UID2

############ if remove the requirement on the number of SM
## calculate number of positive wells defined more than 1 SM
wbcDataFromFile$Pwell = (wbcDataFromFile$SM1>1)+(wbcDataFromFile$SM2>1)
wbcDataFromFile = wbcDataFromFile[wbcDataFromFile$Pwell>=1,]
wbcDataFromFile = wbcDataFromFile[!is.na(wbcDataFromFile$Template),]

## average MAF
wbcDataFromFile$aveMAF = wbcDataFromFile$TotalSM/wbcDataFromFile$TotalUID
maf = c("maf1","maf2")

#wbcDataFromFile = wbcDataFromFile[wbcDataFromFile$TotalUID>=200,]

## use wbc without matched plasma as control, excluding N18s
training = wbcDataFromFile[is.na(wbcDataFromFile$Matched.Plasma), ]
training = training[training$Template!="N18",]
#training = training[!(training$mutID %in% blacklist),]
training$SM3 = training$SM4 = training$SM5 = training$SM6 = NA
training$UID3 = training$UID4 = training$UID5 = training$UID6 = NA
training$maf3 = training$maf4 = training$maf5 = training$maf6 = NA

## normal unmatched WBC mutations
training_wbc.nm = training

## plasma templates with matched WBC
valid = wbcDataFromFile[!is.na(wbcDataFromFile$Matched.Plasma),]
## setting MAF = 0 for wells with less than 100 UIDs
valid$maf1[valid$UID1<100] = 0
valid$maf2[valid$UID2<100] = 0
## unique mutations ID (Plasma Template, mutID)
valid$plsid = paste(valid$Matched.Plasma,valid$mutID)


## calculate maximum MAF in the two wells for each mutation in the plasma
## could use average, instead of max
maxinWBC = ddply(valid, .(plsid), summarise, max = max(maf1,maf2))


#=============
# Plasma Data
#=============
## read data from plasmas
plasmaData = fread(plasmaSeqDataFile)
plasmaData = data.frame(plasmaData)
plasmaData = plasmaData[!(plasmaData$ampMatchName %in% amplicon),]
colnames(plasmaData)[c(25,26,14,16,18,20,22,24)] = c("TotalSM","TotalUID", "UID1","UID2", "UID3","UID4","UID5","UID6")

plasmaData$mutID = paste(plasmaData$Chrom,plasmaData$Position,plasmaData$BaseFrom,plasmaData$BaseTo)

## deleting NA rows and snps
plasmaData = plasmaData[!is.na(plasmaData$Position),]
plasmaData = plasmaData[!(plasmaData$mutID %in% snps),]

## remove duplicated row (same mutation in the same template)
plasmaData$remove = paste(plasmaData$Template,plasmaData$mutID)

dupPairs <- names(which(table(plasmaData$remove)>1))
newuniquenewnormal <- plasmaData%>% 
  filter(remove %in% dupPairs)%>%
  group_by(remove)%>%
  slice(which.max(TotalUID)) %>% 
  as.data.frame()

data_corrected <- rbind(newuniquenewnormal,
                        plasmaData%>% 
                          filter(!(remove %in% dupPairs)))

data_corrected = data.frame(data_corrected)

## numeric columns
col = c("SM1","UID1","SM2","UID2","SM3","UID3","SM4","UID4","SM5","UID5","SM6","UID6","TotalSM","TotalUID")

#colnames(data)[c(17,20,23,26,29,32)] = c("maf1","maf2","maf3","maf4","maf5","maf6")

## change to numeric
data_corrected[,col] = apply(data_corrected[,col], 2, function(x) as.numeric(x))

## calculate MAF
data_corrected$maf1 = data_corrected$SM1/data_corrected$UID1
data_corrected$maf2 = data_corrected$SM2/data_corrected$UID2
data_corrected$maf3 = data_corrected$SM3/data_corrected$UID3
data_corrected$maf4 = data_corrected$SM4/data_corrected$UID4
data_corrected$maf5 = data_corrected$SM5/data_corrected$UID5
data_corrected$maf6 = data_corrected$SM6/data_corrected$UID6

## mark as character
data_corrected$remove = as.character(data_corrected$remove)

## note some templates might only have 2 wells (like normal wbc in this study)
data_corrected$aveMAF = data_corrected$TotalSM/data_corrected$TotalUID
maf = c("maf1","maf2","maf3","maf4","maf5","maf6")
colnames(data_corrected)[6] = "Sample.Category"

## remove water
data_corrected = data_corrected[data_corrected$Template!="Water",]


############# if remove the requirement on the number of SM
## calculating positive wells
data_corrected$Pwell = (data_corrected$SM1>1)+(data_corrected$SM2>1)+(data_corrected$SM3>1)+(data_corrected$SM4>1)+(data_corrected$SM5>1)+(data_corrected$SM6>1)
data_corrected = data_corrected[data_corrected$Pwell>=1,]
data_corrected = data_corrected[!is.na(data_corrected$Template),]

## calculate indel length (SBS is 0)
## not used in this method
data_corrected$indel_length = 0
data_corrected$indel_length[data_corrected$MutType=="Indel" & data_corrected$BaseFrom=="NULL"] = nchar(data_corrected$BaseTo[data_corrected$MutType=="Indel" & data_corrected$BaseFrom=="NULL"]) 
data_corrected$indel_length[data_corrected$MutType=="Indel" & data_corrected$BaseTo=="NULL"] = nchar(data_corrected$BaseFrom[data_corrected$MutType=="Indel" & data_corrected$BaseTo=="NULL"]) 

data_corrected = data_corrected[!(data_corrected$Template %in% c("Water","N18")),]

## all normals 
# LD fix 10/16/2017
nlpls = unique(data_corrected$Template[data_corrected$Sample.Category=="Normal"])

## mutatioins that have matched WBC in the WBC sheets
checkmut = data_corrected[!is.na(data_corrected$WBC.Sample),]
checkmut = checkmut[checkmut$Template %in% valid$Matched.Plasma,]

## setting MAF=0 if less than 100 UIDs
checkmut$maf1[checkmut$UID1<100] = 0
checkmut$maf2[checkmut$UID2<100] = 0
checkmut$maf3[checkmut$UID3<100] = 0
checkmut$maf4[checkmut$UID4<100] = 0
checkmut$maf5[checkmut$UID5<100] = 0
checkmut$maf6[checkmut$UID6<100] = 0

## max MAF in 6 wells from plasmas
maxinPls = ddply(checkmut, .(remove), summarise, max = max(maf1,maf2,maf3,maf4,maf5,maf6))
## corresponding max in matched WBC
maxinPls$WBC = maxinWBC$max[match(maxinPls$remove,maxinWBC$plsid)]

## mutations that are matched in plasma and WBC
maxinPls = maxinPls[!is.na(maxinPls$WBC),]
## calculate the ratio between matched plasma mutation and corresponding WBC mutation
maxinPls$ratio = maxinPls$max/maxinPls$WBC

## "chip" mutations to be removed if the ratio is less than 100
chip_r = 100
## mutations to be removed
chip = maxinPls$remove[maxinPls$ratio<chip_r]


  
  set.seed(1234+r*10)
  ## 10 folds cross validation

  ## cancer type considered in the study
  cancer = c("CRC","Lung","Breast","Pancreas","Ovarian","Esophagus","Liver","Stomach","pancreas")
  
  ## all cancer templates in the cancer category
  #pls.cancer = unique(data_corrected$Template[data_corrected$Category %in% cancer])
  # LD fix 10/16/2017
  pls.cancer = unique(data_corrected$Template[data_corrected$Sample.Category %in% cancer])
  
  # cvRuns-fold cross-validation
  index = cvFolds(length(nlpls), K=cvRuns)
  index2 = cvFolds(length(pls.cancer), K=cvRuns)
  
  # run cross-validation
  for( w in 1:cvRuns)
  {
    ## normals used as control
    nlt = nlpls[index$subsets[index$which!=w]]
    ## cancers used to build the cancer distribution
    plsc.t = pls.cancer[index2$subsets[index2$which!=w]]
    
    ## normal plasma mutations 
    training_nlpls = data_corrected[data_corrected$Template %in% nlt,]
    ## only mutations that are not on blacklist or caution list will be included in the normalization step and the distribution analysis
    
    ## WBC templates
    wbc = unique(training_wbc.nm$Template)
    col = intersect(colnames(training_wbc.nm),colnames(training_nlpls))
    
    ## using the not matched wbc and 80% normal plasma as control to build the distribution for normals
    training = rbind.data.frame(training_wbc.nm[,col],training_nlpls[col])
    
    training$maf1[training$UID1<100] = 0
    training$maf2[training$UID2<100] = 0
    training$maf3[training$UID3<100] = 0
    training$maf4[training$UID4<100] = 0
    training$maf5[training$UID5<100] = 0
    training$maf6[training$UID6<100] = 0
    
    ## eliminate "chip" mutations
    training = training[!(training$remove %in% chip),]
    ## calculate indel length
    training$indel_length = 0
    training$indel_length[training$MutType=="Indel" & training$BaseFrom=="NULL"] = nchar(training$BaseTo[training$MutType=="Indel" & training$BaseFrom=="NULL"]) 
    training$indel_length[training$MutType=="Indel" & training$BaseTo=="NULL"] = nchar(training$BaseFrom[training$MutType=="Indel" & training$BaseTo=="NULL"]) 
    
    ## 1 base pair indel
    indels1 = unique(training$mutID[training$indel_length==1])
    ## more than 1 base pair indel
    indels2 = unique(training$mutID[training$indel_length>=2])
    ## sbs mutations
    sbs = unique(training$mutID[training$indel_length==0])
    
    
    
    ###########  Normalization ###########
    mut = unique(training$mutID)
    ## NB!!! 2 since normal wbc only have 2 wells
    ## create matrix with mutations in normals vs the number of normals times 2 (because of 2 wells)
    cc = 2*length(wbc) + 6*length(nlt)
    mutation_info = matrix(0, nrow = length(mut), ncol = cc)
    
    ## fill in the matrix
    maf_wbc = c("maf1","maf2","maf3","maf4","maf5","maf6")
    for(i in 1:nrow(mutation_info)){
      temp = training[training$mutID==mut[i],]
      temp1 = as.vector(t(temp[,maf_wbc]))
      temp1 = temp1[!is.na(temp1)]
      mutation_info[i,1:(length(temp1))] = temp1
    }
    
    ## count how many times the mutation is present in training. use for normalization only if present more than once
    count = table(c(training$mutID[training$SM1>0], training$mutID[training$SM2>0], training$mutID[training$SM3>0],training$mutID[training$SM4>0],training$mutID[training$SM5>0],training$mutID[training$SM6>0]))
    valid = names(count)[count>1]
    
    ## calculate the mean and max MAF for mutations in normals
    normalization_mutation = data.frame(mutID = mut, stringsAsFactors = FALSE)
    normalization_mutation$mean = 0
    normalization_mutation$max = 0
    for(i in 1:nrow(normalization_mutation)){
      normalization_mutation$mean[i] = mean(mutation_info[i,], na.rm = T)
      normalization_mutation$max[i] = max(mutation_info[i,], na.rm = T)
    }
    
    normalization_mutation = normalization_mutation[normalization_mutation$mean>0,]
    ## normalization is done on all mutations, including blacklist mutations
    
    ## if the mutation has more than 1 positive well in WBCs
    normalization_mutation$valid = (normalization_mutation$mutID %in% valid)
    ## if the mutation is on the blacklist, 1=not on blacklist, 0=yes
    ## normalization_mutation$bl = (!(normalization_mutation$mutID %in% blacklist))
    ## if the mutation should be included for the normalization step, 1=yes, 0=no
    ## normalization_mutation$include = (normalization_mutation$valid)*(normalization_mutation$bl)
    normalization_mutation$include = (normalization_mutation$valid)
    
    normalization_mutation$mutType = ""
    normalization_mutation$mutType[normalization_mutation$mutID %in% indels1] = "Indel1"
    normalization_mutation$mutType[normalization_mutation$mutID %in% indels2] = "Indel2"
    normalization_mutation$mutType[normalization_mutation$mutID %in% sbs] = "SBS"
    
    training$adj1 = training$maf1
    training$adj2 = training$maf2
    training$adj3 = training$maf3
    training$adj4 = training$maf4
    training$adj5 = training$maf5
    training$adj6 = training$maf6
    
    
    ## calculating 25th percentile of the mean MAFs in training mutations (excluding blacklist)
    base = quantile(normalization_mutation$mean[normalization_mutation$include==1],0.25)
    ## not distinguish between indels/SBS
    base1 = base
    base2 = base
    base3 = base
    
    ### distinguish between Indels/SBS (normalized separately)
    #base1 = quantile(normalization_mutation$mean[normalization_mutation$include==1 & normalization_mutation$mutType=="SBS"],0.25)
    #base2 = quantile(normalization_mutation$mean[normalization_mutation$include==1 & normalization_mutation$mutType=="Indel1"],0.25)
    #base3 = quantile(normalization_mutation$mean[normalization_mutation$include==1 & normalization_mutation$mutType=="Indel2"],0.25)
    
    ## calculating adjusted MAF in normals
    for(i in 1:nrow(normalization_mutation)){
      if((normalization_mutation$include[i]==1) & (normalization_mutation$mutType[i]=="SBS")){
        training$adj1[training$mutID==normalization_mutation$mutID[i]] = training$maf1[training$mutID==normalization_mutation$mutID[i]]/(normalization_mutation$mean[i]/base1)
        training$adj2[training$mutID==normalization_mutation$mutID[i]] = training$maf2[training$mutID==normalization_mutation$mutID[i]]/(normalization_mutation$mean[i]/base1)
        training$adj3[training$mutID==normalization_mutation$mutID[i]] = training$maf3[training$mutID==normalization_mutation$mutID[i]]/(normalization_mutation$mean[i]/base1)
        training$adj4[training$mutID==normalization_mutation$mutID[i]] = training$maf4[training$mutID==normalization_mutation$mutID[i]]/(normalization_mutation$mean[i]/base1)
        training$adj5[training$mutID==normalization_mutation$mutID[i]] = training$maf5[training$mutID==normalization_mutation$mutID[i]]/(normalization_mutation$mean[i]/base1)
        training$adj6[training$mutID==normalization_mutation$mutID[i]] = training$maf6[training$mutID==normalization_mutation$mutID[i]]/(normalization_mutation$mean[i]/base1)
      }else if((normalization_mutation$include[i]==1) & (normalization_mutation$mutType[i]=="Indel1")){
        training$adj1[training$mutID==normalization_mutation$mutID[i]] = training$maf1[training$mutID==normalization_mutation$mutID[i]]/(normalization_mutation$mean[i]/base2)
        training$adj2[training$mutID==normalization_mutation$mutID[i]] = training$maf2[training$mutID==normalization_mutation$mutID[i]]/(normalization_mutation$mean[i]/base2)
        training$adj3[training$mutID==normalization_mutation$mutID[i]] = training$maf3[training$mutID==normalization_mutation$mutID[i]]/(normalization_mutation$mean[i]/base2)
        training$adj4[training$mutID==normalization_mutation$mutID[i]] = training$maf4[training$mutID==normalization_mutation$mutID[i]]/(normalization_mutation$mean[i]/base2)
        training$adj5[training$mutID==normalization_mutation$mutID[i]] = training$maf5[training$mutID==normalization_mutation$mutID[i]]/(normalization_mutation$mean[i]/base2)
        training$adj6[training$mutID==normalization_mutation$mutID[i]] = training$maf6[training$mutID==normalization_mutation$mutID[i]]/(normalization_mutation$mean[i]/base2)
      } else if((normalization_mutation$include[i]==1) & (normalization_mutation$mutType[i]=="Indel2")){
        training$adj1[training$mutID==normalization_mutation$mutID[i]] = training$maf1[training$mutID==normalization_mutation$mutID[i]]/(normalization_mutation$mean[i]/base3)
        training$adj2[training$mutID==normalization_mutation$mutID[i]] = training$maf2[training$mutID==normalization_mutation$mutID[i]]/(normalization_mutation$mean[i]/base3)
        training$adj3[training$mutID==normalization_mutation$mutID[i]] = training$maf3[training$mutID==normalization_mutation$mutID[i]]/(normalization_mutation$mean[i]/base3)
        training$adj4[training$mutID==normalization_mutation$mutID[i]] = training$maf4[training$mutID==normalization_mutation$mutID[i]]/(normalization_mutation$mean[i]/base3)
        training$adj5[training$mutID==normalization_mutation$mutID[i]] = training$maf5[training$mutID==normalization_mutation$mutID[i]]/(normalization_mutation$mean[i]/base3)
        training$adj6[training$mutID==normalization_mutation$mutID[i]] = training$maf6[training$mutID==normalization_mutation$mutID[i]]/(normalization_mutation$mean[i]/base3)
      } 
    }
    
    ## build the distribution based on adjusted MAF in controls
    ## no separation between Indel/SBS
    ## here training includes all non-matched WBC and normal plasma controls
    training1 = training
    
    ### Separate between Indel/SBS
    #training1 = training[training$MutType=="SBS",]
    #training2 = training[training$MutType=="Indel" & training$indel_length==1,]
    #training3 = training[training$MutType=="Indel" & training$indel_length>=2,]
    
    data2 = data.frame(perc = unlist(c(training1$adj1,training1$adj2,training1$adj3,training1$adj4,training1$adj5,training1$adj6)))
    data2$UID = unlist(c(training1$UID1,training1$UID2,training1$UID3,training1$UID4,training1$UID5,training1$UID6))
    
    data2 = data2[!is.na(data2$perc),]
    data2$perc = as.numeric(data2$perc)
    
    s1 = 1000
    each = 1000
    
    r1 = data2[data2$UID<s1 & data2$UID>=100,]
    r1$perc.log = as.numeric(log10(r1$perc))
    
    r2 = data2[data2$UID>=s1 & data2$UID<(s1+each*1),]
    r2$perc.log = as.numeric(log10(r2$perc))
    
    r3 = data2[data2$UID>=(s1+each*1) & data2$UID<(s1+each*2),]
    r3$perc.log = as.numeric(log10(r3$perc))
    
    r4 = data2[data2$UID>=(s1+each*2) & data2$UID<(s1+each*3),]
    r4$perc.log = as.numeric(log10(r4$perc))
    
    r5 = data2[data2$UID>=(s1+each*3) & data2$UID<(s1+each*4),]
    r5$perc.log = as.numeric(log10(r5$perc))
    
    r6 = data2[data2$UID>=(s1+each*4) & data2$UID<(s1+each*5),]
    r6$perc.log = as.numeric(log10(r6$perc))
    
    r7 = data2[data2$UID>=(s1+each*5) & data2$UID<(s1+each*6),]
    r7$perc.log = as.numeric(log10(r7$perc))
    
    r8 = data2[data2$UID>=(s1+each*6) & data2$UID<(s1+each*7),]
    r8$perc.log = as.numeric(log10(r8$perc))
    
    r9 = data2[data2$UID>=(s1+each*7) & data2$UID<(s1+each*8),]
    r9$perc.log = as.numeric(log10(r9$perc))
    
    r10 = data2[data2$UID>=(s1+each*8),]
    r10$perc.log = as.numeric(log10(r10$perc))
    
    d1 = density(r1$perc.log)
    d2 = density(r2$perc.log)
    d3 = density(r3$perc.log)
    d4 = density(r4$perc.log)
    d5 = density(r5$perc.log)
    d6 = density(r6$perc.log)
    d7 = density(r7$perc.log)
    d8 = density(r8$perc.log)
    d9 = density(r9$perc.log)
    d10 = density(r10$perc.log)
    
    list_func1 = list(d1, d2, d3, d4, d5, d6, d7,d8,d9,d10)
    list_data1 = list(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10)
    
    l1 = length(list_func1)
    
    
    testSet = data_corrected
    
    testSet = data.frame(testSet)
    ## remove matched plasma mutations and WBC mutations
    testSet = testSet[!(testSet$remove %in% chip),]
    
    ## calculate adjusted MAF for all templates
    testSet$maf1[is.na(testSet$UID1)] = 0
    testSet$maf2[is.na(testSet$UID2)] = 0
    testSet$maf3[is.na(testSet$UID3)] = 0
    testSet$maf4[is.na(testSet$UID4)] = 0
    testSet$maf5[is.na(testSet$UID5)] = 0
    testSet$maf6[is.na(testSet$UID6)] = 0
    
    testSet$maf1[testSet$UID1==0] = 0
    testSet$maf2[testSet$UID2==0] = 0
    testSet$maf3[testSet$UID3==0] = 0
    testSet$maf4[testSet$UID4==0] = 0
    testSet$maf5[testSet$UID5==0] = 0
    testSet$maf6[testSet$UID6==0] = 0
    
    testSet$adj1 = testSet$maf1
    testSet$adj2 = testSet$maf2
    testSet$adj3 = testSet$maf3
    testSet$adj4 = testSet$maf4
    testSet$adj5 = testSet$maf5
    testSet$adj6 = testSet$maf6
    
    testSet$max = 0
    
    ## normalization of mutations in the test set
    for(i in 1:nrow(normalization_mutation)){
      testSet$max[testSet$mutID==normalization_mutation$mutID[i]] = normalization_mutation$max[i]
      if((normalization_mutation$include[i]==1) & (normalization_mutation$mutType[i]=="SBS")){
        testSet$adj1[testSet$mutID==normalization_mutation$mutID[i]] = testSet$maf1[testSet$mutID==normalization_mutation$mutID[i]]/(normalization_mutation$mean[i]/base1)
        testSet$adj2[testSet$mutID==normalization_mutation$mutID[i]] = testSet$maf2[testSet$mutID==normalization_mutation$mutID[i]]/(normalization_mutation$mean[i]/base1)
        testSet$adj3[testSet$mutID==normalization_mutation$mutID[i]] = testSet$maf3[testSet$mutID==normalization_mutation$mutID[i]]/(normalization_mutation$mean[i]/base1)
        testSet$adj4[testSet$mutID==normalization_mutation$mutID[i]] = testSet$maf4[testSet$mutID==normalization_mutation$mutID[i]]/(normalization_mutation$mean[i]/base1)
        testSet$adj5[testSet$mutID==normalization_mutation$mutID[i]] = testSet$maf5[testSet$mutID==normalization_mutation$mutID[i]]/(normalization_mutation$mean[i]/base1)
        testSet$adj6[testSet$mutID==normalization_mutation$mutID[i]] = testSet$maf6[testSet$mutID==normalization_mutation$mutID[i]]/(normalization_mutation$mean[i]/base1)
      }else if((normalization_mutation$include[i]==1) & (normalization_mutation$mutType[i]=="Indel1")){
        testSet$adj1[testSet$mutID==normalization_mutation$mutID[i]] = testSet$maf1[testSet$mutID==normalization_mutation$mutID[i]]/(normalization_mutation$mean[i]/base2)
        testSet$adj2[testSet$mutID==normalization_mutation$mutID[i]] = testSet$maf2[testSet$mutID==normalization_mutation$mutID[i]]/(normalization_mutation$mean[i]/base2)
        testSet$adj3[testSet$mutID==normalization_mutation$mutID[i]] = testSet$maf3[testSet$mutID==normalization_mutation$mutID[i]]/(normalization_mutation$mean[i]/base2)
        testSet$adj4[testSet$mutID==normalization_mutation$mutID[i]] = testSet$maf4[testSet$mutID==normalization_mutation$mutID[i]]/(normalization_mutation$mean[i]/base2)
        testSet$adj5[testSet$mutID==normalization_mutation$mutID[i]] = testSet$maf5[testSet$mutID==normalization_mutation$mutID[i]]/(normalization_mutation$mean[i]/base2)
        testSet$adj6[testSet$mutID==normalization_mutation$mutID[i]] = testSet$maf6[testSet$mutID==normalization_mutation$mutID[i]]/(normalization_mutation$mean[i]/base2)
      }else if((normalization_mutation$include[i]==1) & (normalization_mutation$mutType[i]=="Indel2")){
        testSet$adj1[testSet$mutID==normalization_mutation$mutID[i]] = testSet$maf1[testSet$mutID==normalization_mutation$mutID[i]]/(normalization_mutation$mean[i]/base3)
        testSet$adj2[testSet$mutID==normalization_mutation$mutID[i]] = testSet$maf2[testSet$mutID==normalization_mutation$mutID[i]]/(normalization_mutation$mean[i]/base3)
        testSet$adj3[testSet$mutID==normalization_mutation$mutID[i]] = testSet$maf3[testSet$mutID==normalization_mutation$mutID[i]]/(normalization_mutation$mean[i]/base3)
        testSet$adj4[testSet$mutID==normalization_mutation$mutID[i]] = testSet$maf4[testSet$mutID==normalization_mutation$mutID[i]]/(normalization_mutation$mean[i]/base3)
        testSet$adj5[testSet$mutID==normalization_mutation$mutID[i]] = testSet$maf5[testSet$mutID==normalization_mutation$mutID[i]]/(normalization_mutation$mean[i]/base3)
        testSet$adj6[testSet$mutID==normalization_mutation$mutID[i]] = testSet$maf6[testSet$mutID==normalization_mutation$mutID[i]]/(normalization_mutation$mean[i]/base3)
      }
    }
    
    testSet$perc1 = log10(testSet$adj1)
    testSet$perc2 = log10(testSet$adj2)
    testSet$perc3 = log10(testSet$adj3)
    testSet$perc4 = log10(testSet$adj4)
    testSet$perc5 = log10(testSet$adj5)
    testSet$perc6 = log10(testSet$adj6)
    
    ## cancer templates in cancer distribution
    cancerPT = testSet[testSet$PT.Avg.MAF>5 & (testSet$Template %in% plsc.t),]
    cancerPT = cancerPT[!is.na(cancerPT$Template),]
    cancermaf = data.frame(UID = c(cancerPT$UID1,cancerPT$UID2,cancerPT$UID3,cancerPT$UID4,cancerPT$UID5,cancerPT$UID6), perc = c(cancerPT$adj1,cancerPT$adj2,cancerPT$adj3,cancerPT$adj4,cancerPT$adj5,cancerPT$adj6))
    
    data2 = cancermaf
    data2 = data2[!is.na(data2$perc),]
    data2$perc = as.numeric(data2$perc)
    
    s1 = 1000
    each = 1000
    
    r1 = data2[data2$UID<s1 & data2$UID>=100,]
    r1$perc.log = as.numeric(log10(r1$perc))
    
    r2 = data2[data2$UID>=s1 & data2$UID<(s1+each*1),]
    r2$perc.log = as.numeric(log10(r2$perc))
    
    r3 = data2[data2$UID>=(s1+each*1) & data2$UID<(s1+each*2),]
    r3$perc.log = as.numeric(log10(r3$perc))
    
    r4 = data2[data2$UID>=(s1+each*2) & data2$UID<(s1+each*3),]
    r4$perc.log = as.numeric(log10(r4$perc))
    
    r5 = data2[data2$UID>=(s1+each*3) & data2$UID<(s1+each*4),]
    r5$perc.log = as.numeric(log10(r5$perc))
    
    r6 = data2[data2$UID>=(s1+each*4) & data2$UID<(s1+each*5),]
    r6$perc.log = as.numeric(log10(r6$perc))
    
    r7 = data2[data2$UID>=(s1+each*5) & data2$UID<(s1+each*6),]
    r7$perc.log = as.numeric(log10(r7$perc))
    
    r8 = data2[data2$UID>=(s1+each*6) & data2$UID<(s1+each*7),]
    r8$perc.log = as.numeric(log10(r8$perc))
    
    r9 = data2[data2$UID>=(s1+each*7) & data2$UID<(s1+each*8),]
    r9$perc.log = as.numeric(log10(r9$perc))
    
    r10 = data2[data2$UID>=(s1+each*8),]
    r10$perc.log = as.numeric(log10(r10$perc))
    
    d1 = density(r1$perc.log)
    d2 = density(r2$perc.log)
    d3 = density(r3$perc.log)
    d4 = density(r4$perc.log)
    d5 = density(r5$perc.log)
    d6 = density(r6$perc.log)
    d7 = density(r7$perc.log)
    d8 = density(r8$perc.log)
    d9 = density(r9$perc.log)
    d10 = density(r10$perc.log)
    

    list_func2 = list(d1, d2, d3, d4, d5, d6, d7,d8,d9,d10)
    list_data2 = list(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10)
    
    ## list of data and distribution for normals
    list_func =list(list_func1)
    list_data = list(list_data1)
    
    ## list of data and distribution for cancers
    list_func.c = list(list_func2)
    list_data.c = list(list_data2)
    
    
    testSet$p.value1 = 0.5
    testSet$p.value2 = 0.5
    testSet$p.value3 = 0.5
    testSet$p.value4 = 0.5
    testSet$p.value5 = 0.5
    testSet$p.value6 = 0.5
    
    ## adding the UID group numbers
    testSet$temp1 = apply(testSet[,c("UID1","UID2","UID3","UID4","UID5","UID6")], 1, function(x) min(ceiling((x[1]-(s1-1))/each)+1,l1))
    testSet$temp2 = apply(testSet[,c("UID1","UID2","UID3","UID4","UID5","UID6")], 1, function(x) min(ceiling((x[2]-(s1-1))/each)+1,l1))
    testSet$temp3 = apply(testSet[,c("UID1","UID2","UID3","UID4","UID5","UID6")], 1, function(x) min(ceiling((x[3]-(s1-1))/each)+1,l1))
    testSet$temp4 = apply(testSet[,c("UID1","UID2","UID3","UID4","UID5","UID6")], 1, function(x) min(ceiling((x[4]-(s1-1))/each)+1,l1))
    testSet$temp5 = apply(testSet[,c("UID1","UID2","UID3","UID4","UID5","UID6")], 1, function(x) min(ceiling((x[5]-(s1-1))/each)+1,l1))
    testSet$temp6 = apply(testSet[,c("UID1","UID2","UID3","UID4","UID5","UID6")], 1, function(x) min(ceiling((x[6]-(s1-1))/each)+1,l1))
    
    
    testSet$indel_group = 1
    
    ## calculate the p values
    testSet$adj1_1 = testSet$adj1>0
    testSet$adj2_1 = testSet$adj2>0
    testSet$adj3_1 = testSet$adj3>0
    testSet$adj4_1 = testSet$adj4>0
    testSet$adj5_1 = testSet$adj5>0
    testSet$adj6_1 = testSet$adj6>0
    
    ## calculate p values if the well is from normals
    for(i in 1:l1){
      testSet$p.value1[testSet$temp1==i & testSet$adj1_1==TRUE] = apply(testSet[testSet$temp1==i & testSet$adj1_1==TRUE,c("perc1", "perc2","perc3", "perc4","perc5", "perc6","indel_group")], 1, function(x) return(mean(pnorm(x[1]-list_data[[x[7]]][[i]]$perc.log, sd = list_func[[x[7]]][[i]]$bw, lower.tail=FALSE))))
      testSet$p.value2[testSet$temp2==i & testSet$adj2_1==TRUE] = apply(testSet[testSet$temp2==i & testSet$adj2_1==TRUE,c("perc1", "perc2","perc3", "perc4","perc5", "perc6","indel_group")], 1, function(x) return(mean(pnorm(x[2]-list_data[[x[7]]][[i]]$perc.log, sd = list_func[[x[7]]][[i]]$bw, lower.tail=FALSE))))
      testSet$p.value3[testSet$temp3==i & testSet$adj3_1==TRUE] = apply(testSet[testSet$temp3==i & testSet$adj3_1==TRUE,c("perc1", "perc2","perc3", "perc4","perc5", "perc6","indel_group")], 1, function(x) return(mean(pnorm(x[3]-list_data[[x[7]]][[i]]$perc.log, sd = list_func[[x[7]]][[i]]$bw, lower.tail=FALSE))))
      testSet$p.value4[testSet$temp4==i & testSet$adj4_1==TRUE] = apply(testSet[testSet$temp4==i & testSet$adj4_1==TRUE,c("perc1", "perc2","perc3", "perc4","perc5", "perc6","indel_group")], 1, function(x) return(mean(pnorm(x[4]-list_data[[x[7]]][[i]]$perc.log, sd = list_func[[x[7]]][[i]]$bw, lower.tail=FALSE))))
      testSet$p.value5[testSet$temp5==i & testSet$adj5_1==TRUE] = apply(testSet[testSet$temp5==i & testSet$adj5_1==TRUE,c("perc1", "perc2","perc3", "perc4","perc5", "perc6","indel_group")], 1, function(x) return(mean(pnorm(x[5]-list_data[[x[7]]][[i]]$perc.log, sd = list_func[[x[7]]][[i]]$bw, lower.tail=FALSE))))
      testSet$p.value6[testSet$temp6==i & testSet$adj6_1==TRUE] = apply(testSet[testSet$temp6==i & testSet$adj6_1==TRUE,c("perc1", "perc2","perc3", "perc4","perc5", "perc6","indel_group")], 1, function(x) return(mean(pnorm(x[6]-list_data[[x[7]]][[i]]$perc.log, sd = list_func[[x[7]]][[i]]$bw, lower.tail=FALSE))))
    }
    
    ## setting p values to 1-10^-8 if 0 SM
    testSet$p.value1[testSet$adj1_1==FALSE] = 1-10^-8
    testSet$p.value2[testSet$adj2_1==FALSE] = 1-10^-8
    testSet$p.value3[testSet$adj3_1==FALSE] = 1-10^-8
    testSet$p.value4[testSet$adj4_1==FALSE] = 1-10^-8
    testSet$p.value5[testSet$adj5_1==FALSE] = 1-10^-8
    testSet$p.value6[testSet$adj6_1==FALSE] = 1-10^-8
    
    ## pvalues for wells with less than 100 UIDs are calculate using binomial, p is the aveMAF from 6 wells
    testSet$p.value1[!is.na(testSet$UID1) & testSet$UID1<100] = apply(testSet[!is.na(testSet$UID1) & testSet$UID1<100, c("SM1","maf1","UID1","aveMAF")], 1, function(x) return(pbinom(max(0,x[1]-1), x[3], x[4],lower.tail = FALSE)))
    testSet$p.value2[!is.na(testSet$UID2) & testSet$UID2<100] = apply(testSet[!is.na(testSet$UID2) & testSet$UID2<100, c("SM2","maf2","UID2","aveMAF")], 1, function(x) return(pbinom(max(0,x[1]-1), x[3], x[4],lower.tail = FALSE)))
    testSet$p.value3[!is.na(testSet$UID3) & testSet$UID3<100] = apply(testSet[!is.na(testSet$UID3) & testSet$UID3<100, c("SM3","maf3","UID3","aveMAF")], 1, function(x) return(pbinom(max(0,x[1]-1), x[3], x[4],lower.tail = FALSE)))
    testSet$p.value4[!is.na(testSet$UID4) & testSet$UID4<100] = apply(testSet[!is.na(testSet$UID4) & testSet$UID4<100, c("SM4","maf4","UID4","aveMAF")], 1, function(x) return(pbinom(max(0,x[1]-1), x[3], x[4],lower.tail = FALSE)))
    testSet$p.value5[!is.na(testSet$UID5) & testSet$UID5<100] = apply(testSet[!is.na(testSet$UID5) & testSet$UID5<100, c("SM5","maf5","UID5","aveMAF")], 1, function(x) return(pbinom(max(0,x[1]-1), x[3], x[4],lower.tail = FALSE)))
    testSet$p.value6[!is.na(testSet$UID6) & testSet$UID6<100] = apply(testSet[!is.na(testSet$UID6) & testSet$UID6<100, c("SM6","maf6","UID6","aveMAF")], 1, function(x) return(pbinom(max(0,x[1]-1), x[3], x[4],lower.tail = FALSE)))
    
    ## 0 wells with <100 UIDs are assigned p values 1-10^-8 (theoretically it is 1)
    testSet$p.value1[!is.na(testSet$UID1) & testSet$UID1<100 & testSet$SM1==0] = 1-10^-8
    testSet$p.value2[!is.na(testSet$UID2) & testSet$UID2<100 & testSet$SM2==0] = 1-10^-8
    testSet$p.value3[!is.na(testSet$UID3) & testSet$UID3<100 & testSet$SM3==0] = 1-10^-8
    testSet$p.value4[!is.na(testSet$UID4) & testSet$UID4<100 & testSet$SM4==0] = 1-10^-8
    testSet$p.value5[!is.na(testSet$UID5) & testSet$UID5<100 & testSet$SM5==0] = 1-10^-8
    testSet$p.value6[!is.na(testSet$UID6) & testSet$UID6<100 & testSet$SM6==0] = 1-10^-8
    
    ## p values for wells with 0 UID are assigned as 0.5 (when inverting it has 0)
    testSet$p.value1[testSet$UID1==0] = 0.5
    testSet$p.value2[testSet$UID2==0] = 0.5
    testSet$p.value3[testSet$UID3==0] = 0.5
    testSet$p.value4[testSet$UID4==0] = 0.5
    testSet$p.value5[testSet$UID5==0] = 0.5
    testSet$p.value6[testSet$UID6==0] = 0.5
    
    
    ### calculating p values if the well is from cancers
    testSet$p.value1.c = 0.5
    testSet$p.value2.c = 0.5
    testSet$p.value3.c = 0.5
    testSet$p.value4.c = 0.5
    testSet$p.value5.c = 0.5
    testSet$p.value6.c = 0.5
    
    
    for(i in 1:l1){
      testSet$p.value1.c[testSet$temp1==i & testSet$adj1_1==TRUE] = apply(testSet[testSet$temp1==i & testSet$adj1_1==TRUE,c("perc1", "perc2","perc3", "perc4","perc5", "perc6","indel_group")], 1, function(x) return(mean(pnorm(x[1]-list_data.c[[x[7]]][[i]]$perc.log, sd = list_func.c[[x[7]]][[i]]$bw, lower.tail=FALSE))))
      testSet$p.value2.c[testSet$temp2==i & testSet$adj2_1==TRUE] = apply(testSet[testSet$temp2==i & testSet$adj2_1==TRUE,c("perc1", "perc2","perc3", "perc4","perc5", "perc6","indel_group")], 1, function(x) return(mean(pnorm(x[2]-list_data.c[[x[7]]][[i]]$perc.log, sd = list_func.c[[x[7]]][[i]]$bw, lower.tail=FALSE))))
      testSet$p.value3.c[testSet$temp3==i & testSet$adj3_1==TRUE] = apply(testSet[testSet$temp3==i & testSet$adj3_1==TRUE,c("perc1", "perc2","perc3", "perc4","perc5", "perc6","indel_group")], 1, function(x) return(mean(pnorm(x[3]-list_data.c[[x[7]]][[i]]$perc.log, sd = list_func.c[[x[7]]][[i]]$bw, lower.tail=FALSE))))
      testSet$p.value4.c[testSet$temp4==i & testSet$adj4_1==TRUE] = apply(testSet[testSet$temp4==i & testSet$adj4_1==TRUE,c("perc1", "perc2","perc3", "perc4","perc5", "perc6","indel_group")], 1, function(x) return(mean(pnorm(x[4]-list_data.c[[x[7]]][[i]]$perc.log, sd = list_func.c[[x[7]]][[i]]$bw, lower.tail=FALSE))))
      testSet$p.value5.c[testSet$temp5==i & testSet$adj5_1==TRUE] = apply(testSet[testSet$temp5==i & testSet$adj5_1==TRUE,c("perc1", "perc2","perc3", "perc4","perc5", "perc6","indel_group")], 1, function(x) return(mean(pnorm(x[5]-list_data.c[[x[7]]][[i]]$perc.log, sd = list_func.c[[x[7]]][[i]]$bw, lower.tail=FALSE))))
      testSet$p.value6.c[testSet$temp6==i & testSet$adj6_1==TRUE] = apply(testSet[testSet$temp6==i & testSet$adj6_1==TRUE,c("perc1", "perc2","perc3", "perc4","perc5", "perc6","indel_group")], 1, function(x) return(mean(pnorm(x[6]-list_data.c[[x[7]]][[i]]$perc.log, sd = list_func.c[[x[7]]][[i]]$bw, lower.tail=FALSE))))
    }
    
    
    testSet$p.value1.c[testSet$adj1_1==FALSE] = 1-10^-8
    testSet$p.value2.c[testSet$adj2_1==FALSE] = 1-10^-8
    testSet$p.value3.c[testSet$adj3_1==FALSE] = 1-10^-8
    testSet$p.value4.c[testSet$adj4_1==FALSE] = 1-10^-8
    testSet$p.value5.c[testSet$adj5_1==FALSE] = 1-10^-8
    testSet$p.value6.c[testSet$adj6_1==FALSE] = 1-10^-8
    
    ## pvalues for wells with less than 100 UIDs are calculate using binomial
    testSet$p.value1.c[!is.na(testSet$UID1) & testSet$UID1<100] = apply(testSet[!is.na(testSet$UID1) & testSet$UID1<100, c("SM1","maf1","UID1","aveMAF")], 1, function(x) return(pbinom(max(0,x[1]-1), x[3], x[4],lower.tail = FALSE)))
    testSet$p.value2.c[!is.na(testSet$UID2) & testSet$UID2<100] = apply(testSet[!is.na(testSet$UID2) & testSet$UID2<100, c("SM2","maf2","UID2","aveMAF")], 1, function(x) return(pbinom(max(0,x[1]-1), x[3], x[4],lower.tail = FALSE)))
    testSet$p.value3.c[!is.na(testSet$UID3) & testSet$UID3<100] = apply(testSet[!is.na(testSet$UID3) & testSet$UID3<100, c("SM3","maf3","UID3","aveMAF")], 1, function(x) return(pbinom(max(0,x[1]-1), x[3], x[4],lower.tail = FALSE)))
    testSet$p.value4.c[!is.na(testSet$UID4) & testSet$UID4<100] = apply(testSet[!is.na(testSet$UID4) & testSet$UID4<100, c("SM4","maf4","UID4","aveMAF")], 1, function(x) return(pbinom(max(0,x[1]-1), x[3], x[4],lower.tail = FALSE)))
    testSet$p.value5.c[!is.na(testSet$UID5) & testSet$UID5<100] = apply(testSet[!is.na(testSet$UID5) & testSet$UID5<100, c("SM5","maf5","UID5","aveMAF")], 1, function(x) return(pbinom(max(0,x[1]-1), x[3], x[4],lower.tail = FALSE)))
    testSet$p.value6.c[!is.na(testSet$UID6) & testSet$UID6<100] = apply(testSet[!is.na(testSet$UID6) & testSet$UID6<100, c("SM6","maf6","UID6","aveMAF")], 1, function(x) return(pbinom(max(0,x[1]-1), x[3], x[4],lower.tail = FALSE)))
    
    ## 0 wells with <100 UIDs are assigned p values 1-10^-8 (theoretically it is 1)
    testSet$p.value1.c[!is.na(testSet$UID1) & testSet$UID1<100 & testSet$SM1==0] = 1-10^-8
    testSet$p.value2.c[!is.na(testSet$UID2) & testSet$UID2<100 & testSet$SM2==0] = 1-10^-8
    testSet$p.value3.c[!is.na(testSet$UID3) & testSet$UID3<100 & testSet$SM3==0] = 1-10^-8
    testSet$p.value4.c[!is.na(testSet$UID4) & testSet$UID4<100 & testSet$SM4==0] = 1-10^-8
    testSet$p.value5.c[!is.na(testSet$UID5) & testSet$UID5<100 & testSet$SM5==0] = 1-10^-8
    testSet$p.value6.c[!is.na(testSet$UID6) & testSet$UID6<100 & testSet$SM6==0] = 1-10^-8
    
    ## p values for wells with 0 UID are assigned as 0.5 (when inverting it is 0)
    testSet$p.value1.c[testSet$UID1==0] = 0.5
    testSet$p.value2.c[testSet$UID2==0] = 0.5
    testSet$p.value3.c[testSet$UID3==0] = 0.5
    testSet$p.value4.c[testSet$UID4==0] = 0.5
    testSet$p.value5.c[testSet$UID5==0] = 0.5
    testSet$p.value6.c[testSet$UID6==0] = 0.5
    
    ## ratio p.cancer/p.normal
    testSet$r1 = testSet$p.value1.c/testSet$p.value1
    testSet$r2 = testSet$p.value2.c/testSet$p.value2
    testSet$r3 = testSet$p.value3.c/testSet$p.value3
    testSet$r4 = testSet$p.value4.c/testSet$p.value4
    testSet$r5 = testSet$p.value5.c/testSet$p.value5
    testSet$r6 = testSet$p.value6.c/testSet$p.value6
    
    
    testSet$blacklist = testSet$mutID %in% blacklist	
    
    result = testSet[testSet$TotalUID>=200,]
    result = result[!(result$Template %in% c(nlt,"Water","N18",plsc.t)),]
    result$Template = as.character(result$Template)
    
    save.image(file = paste0("PvalueRatio1209cv_5perc_bl_",r,"_", w,"_1006data.rda"))
  }
  




