# tissue recognition analysis
# protein data from 12-09-17
library(xlsx)
library(randomForest)
library(bmrm)
library(caret)
library(dplyr)
library(data.table)
library(Boruta)

#=============	
# runs RF for tissue recognition with 10-fold CV.
# use predefined iterations and training-test assignments
# Tissue is in sampleType column of the dat object
# keep cancers with more than 150 as is, for those that are less, repeat them several times and sample the rest
runRF_balancedSample = function(dat, fold, BalanceSampleCount = 150)
{
	require(bmrm)
	require(randomForest)
	require(caret)
	require(dplyr)
#	NIter = unique(iterFold$iteration)
#print(head(dat))
	AllpredictionsNIter_tissue = c()
	  MyIndices <- fold
	  nTissue = nlevels(as.factor(dat$sampleType))
	  Allpredictions = vector(length = length(MyIndices)) 
	  Allpredictions_votes = matrix(nrow = length(MyIndices),ncol = nTissue)
	  
	  for( i in unique(MyIndices)){
		trainSamples <- which(MyIndices!=i)
		testSamples <- which(MyIndices==i)
#print(table(dat$sampleType[trainSamples]))
    trainSamplesResampled <- tapply(trainSamples,
                    dat$sampleType[trainSamples],
                    FUN = function(x)
                      if (length(x)>BalanceSampleCount)x else {
						s = rep(x, floor(BalanceSampleCount/length(x)));
						return(c(s,sample(x,size = BalanceSampleCount-length(s))))})  %>%
                  unlist()
    
    ThisRF <- randomForest(as.factor(sampleType)~.,
                           data =dat[trainSamplesResampled,],
						   cutoff = rep(1/nTissue,nTissue))

	Allpredictions[testSamples] <- predict(ThisRF,dat[testSamples,]) %>% as.character()

	Allpredictions_votes[testSamples,] <- predict(ThisRF,dat[testSamples,], type = 'vote')
		
	  }
	  names(Allpredictions) = rownames(Allpredictions_votes)=rownames(dat)
	  colnames(Allpredictions_votes) = levels(as.factor(dat$sampleType))
	  conf = confusionMatrix(Allpredictions,reference = dat$sampleType)
	  AllpredictionsNIter_tissue <- list(prediction = Allpredictions, 
		confusion = conf, votes = Allpredictions_votes)
	print(conf$overall[1])
	return(AllpredictionsNIter_tissue)
}

# extracts positive samples and folds for an iteration and runs RF for this iteration
runOneIter = function(obj, colNames = colnames(obj), addInfo = NULL, tumSamp, sampType)
{
# positive samples only
	bahmanSamp = rownames(obj)[which(obj$LogitCall)] # 
	s = intersect(tumSamp,rownames(obj)[which(obj$LogitCall)]) # 
	print(length(s))
	if (!is.null(addInfo)) prot = data.frame(sampleType = sampType[s], obj[s,colNames],addInfo[s,]) else 
		prot = cbind(sampleType = obj[s,'Tissue'],obj[s,colNames])
	return(runRF_balancedSample(prot, fold = obj[s,'fold']))
}

#===================================================
# loading results from cancer detection
load('Results/ForLuda8proteins.rda')

#==============================
# read in protein data
protData = read.xlsx2(file = 'MEGA Protein Dataset for Random Forest, December 9, 2017.xlsx', sheetIndex = 1, row.names = 1, check.names = F, stringsAsFactors = F)

# remove Vitronectin, LRG-1, Stage, Sample Type
protOnly = protData %>% select(-`Vitronectin`, -`LRG-1`, -`Stage`, -`Sample Type`)

protOnly = sapply(protOnly, as.numeric)
#protOnly = t(do.call('rbind',protOnly))
rownames(protOnly) = rownames(protData)
colnames(protOnly) = gsub('-','',colnames(protOnly))
colnames(protOnly) = gsub(' ','',colnames(protOnly))
protOnly[is.na(protOnly)] = 0

#==================================
# create tumor type object
sampleType = as.character(protData$`Sample Type`)
# combine esophagus and stomach cancers to one group
sampleType[which(sampleType%in%c('Esophagus','Stomach'))] = 'st_eso'
names(sampleType) = rownames(protData)

# IDs of tumor and normal samples
tumSamp = names(sampleType)[which(sampleType != 'Normal')] # 1005
normSamp = names(sampleType)[which(sampleType == 'Normal')] # 812

#====================================
# patient info
sampInfo <- read.xlsx2("aar3247_Suppl. Excel_seq3_v1.xlsx",sheetIndex = 4,
                                           stringsAsFactors = F, startRow = 3, endRow = 1820)
# gender
gender = sampInfo[c("Plasma.sample.ID..","Sex")]
gender = gender[!duplicated(gender[,1]),]
n = gender[,1]
gender = gender[,2]
names(gender) = n
gender[which(gender %in% c('Female','female'))] = 'F'
gender[which(gender %in% c('Male','male'))] = 'M'
gender[which(gender %in% c('unknown',''))] = NA

#=============================
# add matrix with omegas per gene per patient.
#=============================
# make 10 matrices from Analysis_1
omegaPerPatient = omegaPerPatientK13 = vector(mode = 'list',length = 10)
s = c(tumSamp, normSamp)
for(i in 1:10)
{
	load(paste0('allResults_fromLuMethod_',i,'_20171209.rda'))
	perPat = split(allResults,allResults$Template)
	perPat = lapply(perPat, function(x) x[which(!x[,'blacklist']),c('Gene','omega')])
	# get maximum omega for each gene
	perPatMaxOmega = sapply(perPat, function(x) tapply(x$omega, INDEX = x$Gene, FUN = max, na.rm = T))
	genes = unique(allResults$Gene)
	s1 = union(s,names(perPatMaxOmega))
	omegaPerPatient[[i]] = matrix(0,ncol = length(genes), nrow = length(s1), dimnames = list(s1,genes))
	for(j in names(perPatMaxOmega)) omegaPerPatient[[i]][j,names(perPatMaxOmega[[j]])] = perPatMaxOmega[[j]]
	
	# separate KRAS codon 13
	resultKRAS13 = allResults %>% filter(Gene == 'KRAS', codonStart == 13)
	resultNoK13 = allResults %>% filter(Gene != 'KRAS' | codonStart != 13)
	resultKRAS13[,'Gene'] = 'KRAS13'
	allResults = rbind(resultNoK13,resultKRAS13)
	# get maximum omega for each gene
	perPatMaxOmega = sapply(perPat, function(x) tapply(x$omega, INDEX = x$Gene, FUN = max, na.rm = T))
	genes = unique(allResults$Gene)
	s1 = union(s,names(perPatMaxOmega))
	omegaPerPatientK13[[i]] = matrix(0,ncol = length(genes), nrow = length(s1), dimnames = list(s1,genes))
	for(j in names(perPatMaxOmega)) omegaPerPatientK13[[i]][j,names(perPatMaxOmega[[j]])] = perPatMaxOmega[[j]]
}

#===========================
# run all proteins, gender and mutations
set.seed(12485)
res8 = vector(length = 10)
for(i in 1:10)
{
	res8[i] = list(runOneIter(cbind(protOnly,ForLuda[[i]][rownames(protOnly),c('LogitCall','fold')]),
			colNames = 1:39, 
			addInfo = data.frame(gender = gender[rownames(protOnly)],omegaPerPatient[[i]][rownames(protOnly),]), 
			tumSamp = tumSamp, sampType = sampleType))
}

# save results for the paper
tab = sapply(res8, FUN = function(x) c(x$confusion$overall[1], x$confusion$byClass[,1])) 
tab = cbind(tab,mean = apply(tab,1,mean))
# mean accuracy 0.6968313
xlsFile = paste('Results/tissueRecog_results_',Sys.Date(),'.xlsx')
write.xlsx(tab, file = xlsFile, sheetName = 'prot +  mut')
# confusion matrix with fractions
conf = res8[[4]]$confusion$table
write.csv(sweep(conf, 2, apply(conf,2,sum),'/'), file = paste('Results/tissueRecog_results_',Sys.Date(),'.csv'), quote = F)
# per samples votes
tab = cbind(res8[[4]]$votes, prediction = res8[[4]]$prediction, tissue = sampleType[names(res8[[4]]$prediction)])
write.xlsx(tab, file = xlsFile, sheetName = 'Iteration 4. Votes', append = T)

