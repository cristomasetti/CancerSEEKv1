library(tibble)
library(dplyr)
library(bmrm)
library(pROC)
library(glmnet)
library(xlsx)
library(dplyr)
#Functions

writeCSVMultiple <- function(x,file){
  if(dir.exists("file")!=F){
    cat("Creating folder",file,"\n")
    dir.create(file)
  }
  for(i in seq_along(x)){
    write.csv(x = x[[i]],
              file = paste(file, "_",names(x)[[i]],".csv",sep = ""))
  }
}


QFunclog <- function(x,Threshold) 
  ifelse(x>Threshold,log10(x+1),0)
QFuncjumpramp <- function(x,Threshold) 
  ifelse(x>Threshold,x,0)


printRemoveDropResults <- function(ThisResult,whichRound="Round 4"){
  cat("Sensitivity:",(mean(unlist(ThisResult$AllSensitivity))),"\n")
  
  print("Only this protein:")
  print((rownames(ThisResult$AllProbPredictionsMinusFeatureNIter) %>% 
           sapply(FUN = function(x) ThisResult$AllProbPredictionsMinusFeatureNIter[[x,"LR"]]%>% 
                    dplyr::select(-BV) %>%
                    summarise_all(.funs = funs(mean(.[which(CancerStatus==T)]>ThisResult$AllThresh99NIter[[x,"LR"]]))))) %>% 
          apply(MARGIN = c(1,2),as.numeric) %>% 
          rowMeans()%>% sort()%>% t()%>% t())
  
  
  print("Remove this protein:")
  print((rownames(ThisResult$AllProbPredictionsMinusFeatureNIter) %>% 
           sapply(FUN = function(x) ThisResult$AllProbPredictionsDropFeatureNIter[[x,"LR"]]%>% 
                    dplyr::select(-BV) %>%
                    summarise_all(.funs = funs(mean(.[which(CancerStatus==T)]>ThisResult$AllThresh99NIter[[x,"LR"]]))))) %>% 
          apply(MARGIN = c(1,2),as.numeric) %>% 
          rowMeans()%>%sort()%>%t()%>%t())
  
  
  mycoefs <- ThisResult$Alllogit %>% coef(s=0) %>% abs() 
  FirstIteration <- ThisResult$Master4ClassificationwithClassAll[whichRound,"LR"][[1]] 
  MinDiff <- FirstIteration  %>% 
    mutate_at(.vars = setdiff(colnames(FirstIteration),c("CancerStatus","maxOmega","fold")),
              .funs=funs(QFuncjumpramp(.,quantile(.[which(CancerStatus==T)],probs=0.95)))) %>% 
    summarise_at(.vars = setdiff(colnames(FirstIteration),
                                 c("CancerStatus","fold")),
                 .funs = funs(mean(.[which(CancerStatus==T)])-mean(.[which(CancerStatus==F)]))) %>%
    as.data.frame() 
  
  
  RankingCoefs <- data.frame('Ranking Score'= as.vector(MinDiff*mycoefs[names(MinDiff),1]) %>%
                               t(),
                             'Coefs' = mycoefs[names(MinDiff),1] %>% 
                               t()
                             %>%t())
  
  print(RankingCoefs%>% 
          rownames_to_column("Protein")%>%
          arrange(Ranking.Score))
  
  
}

ProbsEachRound <-  function(ThisResult,csvFile=NULL){
  Probs = ThisResult$AllProbPredictionsNIter[,"LR"] %>% sapply(FUN = function(x) x)
  Probs <- rbind(Probs,Thresh=ThisResult$AllThresh99NIter[,"LR"] %>% sapply(FUN = function(x) x))
  Probs <- cbind(Probs,
                 CancerStatus =
                   ThisResult$Master4ClassificationwithClassAll[["Round 1","LR"]][rownames(Probs),"CancerStatus"])
  
  if(!is.null(csvFile))
    write.csv(file = csvFile,x = Probs)
  return(Probs)
}

Generate4LudaObject <- function(ThisResult,rdaFile,TissueType=TissueType_RFData){
  ForLuda <- sapply(setNames(paste("Round",1:10),paste("Round",1:10)),
                     function(x) cbind(ThisResult$Master4ClassificationwithClassAll[[x,"LR"]],
                                       LogitCall=ThisResult$AllpredictionsNIter[[x,"LR"]],
                                       LogitProb=ThisResult$AllProbPredictionsNIter[[x,"LR"]],
                                       Tissue=TissueType_RFData[rownames(ThisResult$Master4ClassificationwithClassAll[[x,"LR"]])
                                         
                                       ]),simplify = F)
  ForLudaFeatureRanks <- ThisResult$AlllogitFeaturesRank
  ForLudaSensitivity <- ThisResult$AllSensitivity
  
  
  save(file=rdaFile,
       list= c("ForLuda","ForLudaFeatureRanks","ForLudaSensitivity","ThisResult"))
  
  return( list(ForLuda=ForLuda,
             ForLudaFeatureRanks=ForLudaFeatureRanks,
             ForLudaSensitivity=ForLudaSensitivity,
             ThisResult=ThisResult))
  
  
}


cv.CancervsNormal.QuantileReplace.Mutatation.knownCV <- function(CancerStatus,
                                                                 # a vector containing cancer status (T=Cancer)
                                                                 Proteins,
                                                                 # a matrix of protien vector
                                                                 Mutations,
                                                                 # a matrix of Mutations
                                                                 # It must at include BV.,
                                                                 # iteration , fold
                                                                 
                                                                 #Methods=c("RF","Logit"),
                                                                 # Methods to test
                                                                 ProteinsHigherCancer = T,
                                                                 #If true, we only keep proteins with higher median in cancer in the sense of wilcoxon
                                                                 ForcePositiveCoeff = F,
                                                                 #If true, we force the coefficients to be positive (This is togher measure than ProteinsHigherCancer)
                                                                 #myseed=1234,
                                                                 #NIter=10,
                                                                 #nfold=10,
                                                                 SpecSet = 0.99, #Specificity we shoot for
                                                                 FeaturesToQuantile, # Features to quantile, set to Null if you do not wish quantiling
                                                                 QuantileF=0.95,
                                                                 QuantileFunction= function(x,Threshold) 
                                                                   ifelse(x>Threshold,x,0),
                                                                 NoPenaltyonCtdna=F,
                                                                 NoLassoPenaly = F){
  
  
  
  mylower <- ifelse(ForcePositiveCoeff,0,-Inf)
  if(ProteinsHigherCancer){
    #apply Wilcoxon to filter proteins if ProteinsHigherCancer==T
    
    Proteins2Use <- (apply(Proteins,MARGIN = 2,
                           FUN = function(x) wilcox.test(
                             x[which(CancerStatus==T)],
                             x[which(CancerStatus==F)],
                             alternative = "greater")$p.value)<0.05) %>%
      which() %>% 
      names()
    
    
    # Proteins2Use <- (colMeans(Proteins[which(CancerStatus==T),,drop=F])- 
    #                    colMeans(Proteins[which(CancerStatus==F),,drop=F])>0) %>% 
    #   which() %>% 
    #   names()
    
  }else{
    Proteins2Use <- colnames(Proteins)
  }
  
  
  # if(is.null(Mutations)){
  #   Master4ClassificationwithClass <- data.frame(CancerStatus=CancerStatus,
  #                                                Proteins[,Proteins2Use,drop=F])
  # }else{
  
  
  #}
  
  #Remove the penalty for ctdna if NoPenaltyonCtdna==T
  if(NoPenaltyonCtdna){
    PenaltyFactor <- c(rep(1,length(Proteins2Use)),rep(0,ncol(Mutations)-3))
  }else{
    PenaltyFactor <-c(rep(1,length(Proteins2Use)),rep(1,ncol(Mutations)-3))
  }
  
  #Only FeaturesToQuantile which have survived the Wilcoxon are being quantiles
  if(!is.null(FeaturesToQuantile))
    FeaturesToQuantile <- intersect(FeaturesToQuantile,Proteins2Use)
  
  #RFFeatures <- colnames(Master4ClassificationwithClass) %>% setdiff("CancerStatus")
  
  
  #augmenting the quantiled and non-quantiled proteins
  MutationsRound1st <- Mutations %>% 
    filter(iteration==1) %>%
    column_to_rownames("BV") 
  MutationsRound1st <- MutationsRound1st[rownames(Proteins),]
  MutationsRound1st[is.na(MutationsRound1st)] <- 0
  
  Master4ClassificationwithClass <- data.frame(CancerStatus=CancerStatus,
                                               Proteins[,Proteins2Use,drop=F],
                                               MutationsRound1st %>% 
                                                 dplyr::select(-iteration,-fold))
  if(is.null(FeaturesToQuantile)){
    Master4ClassificationwithClassAppQuantAug <-
      Master4ClassificationwithClass
  }else{
    
    
    ThrehsoldsTraining <- Master4ClassificationwithClass %>% 
      filter(CancerStatus==F) %>%
      summarise_at(.vars = FeaturesToQuantile,# to implement ramp you can manipulate here
                   .funs = funs(quantile(.,probs = QuantileF))) 
    
    Master4ClassificationwithClassAppQuantAug <- # to implement ramp you can manipulate here
      Master4ClassificationwithClass %>% 
      mutate_(.dots = setNames(paste("QuantileFunction(",
                                     names(ThrehsoldsTraining),
                                     ",",
                                     ThrehsoldsTraining,")",
                                     sep = ""),
                               names(ThrehsoldsTraining)))
    
  } 
  
  #Do cross-validation without quantiling for the Logit
  Alllogit <- cv.glmnet(x= Master4ClassificationwithClassAppQuantAug %>% 
                          dplyr::select(-CancerStatus) %>%
                          as.matrix(),
                        y=Master4ClassificationwithClass$CancerStatus,
                        type.measure = "auc",
                        family="binomial",
                        nfolds = length(unique(MutationsRound1st$fold)),
                        penalty.factor = PenaltyFactor,
                        lower.limits=mylower,
                        keep=TRUE)
  AllLogitThreshod <- apply(Alllogit$fit.preval[which(Master4ClassificationwithClass$CancerStatus==F),],
                            MARGIN = 2,
                            FUN = quantile,probs=SpecSet,na.rm=T) %>% 
    as.numeric() ## to implement ramp you can manipulate here
  
  AllSensity <- sapply(seq_along(AllLogitThreshod),
                       FUN = function(x) 
                         mean(Alllogit$fit.preval[ which(Master4ClassificationwithClass$CancerStatus==T),x]
                              >AllLogitThreshod[x]))
  #find the best lambda suited for desired SpecSet
  if(NoLassoPenaly){
    lambda.best <- 0
  }else{
    lambda.best <-  Alllogit$lambda[which.max(AllSensity)]
  }
  AlllogitFeatures <- which(as.matrix(coef(Alllogit,s=lambda.best))!=0,
                            arr.ind = T,
                            useNames = T) %>% 
    rownames() %>% 
    setdiff("(Intercept)")
  
  NIter <- unique(Mutations$iteration)
  
  AllProbPredictionsNIter <- matrix(data = list(),
                                    ncol = 1,
                                    nrow = length(NIter),
                                    dimnames = list(paste("Round",NIter),
                                                    "LR"))
  AllThresh99NIter <- AllProbPredictionsNIter
  
  AllpredictionsNIter <- AllProbPredictionsNIter
  AllProbPredictionsNIter <- AllProbPredictionsNIter
  AllSensitivity <- AllProbPredictionsNIter
  Master4ClassificationwithClassAll <- AllProbPredictionsNIter
  AllProbPredictionsMinusFeatureNIter <- AllProbPredictionsNIter
  AllProbPredictionsDropFeatureNIter <- AllProbPredictionsNIter
  #set.seed(myseed)
  #AllMyIndices <- sapply(1:NIter,FUN = function(x)
  #  balanced.cv.fold(Master4ClassificationwithClass$CancerStatus) %>% as.numeric,simplify = F) 
  
  #names(AllMyIndices)<- paste("Round",NIter,sep = "")
  
  set.seed(1234)
  for( j in NIter){
    #MyIndices <- AllMyIndices[[j]]
    ThisRound <-  Mutations %>%
      filter(iteration==j) 
    
    uniqueMyIndices <- ThisRound %>% 
      dplyr::select(fold) %>% 
      unlist() %>% 
      unique()
    
    Allpredictions <- setNames(vector(length = nrow(Proteins)),rownames(Proteins))
    ProbPredictions <- setNames(vector(length = nrow(Proteins),mode = "numeric") ,rownames(Proteins))
    ProbPredictionsMinusFeature <- matrix(0,
                                          nrow=nrow(Proteins),
                                          ncol= ncol(Master4ClassificationwithClass)-1 ,
                                          dimnames = list(rownames(Proteins),Master4ClassificationwithClass %>% colnames() %>% setdiff("CancerStatus"))
                                          
                                          
    )
    ProbPredictionsDropFeature <- ProbPredictionsMinusFeature
    #for(k in seq_along(Methods)){
    MutationsThisRound <- ThisRound %>%
      column_to_rownames("BV")
    MutationsThisRound <- MutationsThisRound[rownames(Proteins),]
    MutationsThisRound[is.na(MutationsThisRound)] <- 0
    
    Master4ClassificationwithClass <- data.frame(CancerStatus=CancerStatus,
                                                 Proteins[,Proteins2Use,drop=F],
                                                 MutationsThisRound %>% 
                                                   dplyr::select(-iteration,-fold))
    for( i in uniqueMyIndices){
      trainSamples <-  ThisRound %>% #which(MyIndices!=i)
        filter(fold!=i) %>%
        dplyr::select(BV) %>%
        unlist()
      testSamples <- ThisRound %>% #which(MyIndices!=i)
        filter(fold==i) %>%
        dplyr::select(BV) %>%
        unlist()
      
      
      if(is.null(FeaturesToQuantile)){
        Master4ClassificationwithClassAppQuantAug <-
          Master4ClassificationwithClass
      }else{
        
        ThrehsoldsTraining <- Master4ClassificationwithClass[trainSamples,] %>% 
          filter(CancerStatus==F) %>%
          summarise_at(.vars = FeaturesToQuantile,# to implement ramp you can manipulate here
                       .funs = funs(quantile(.,probs = QuantileF))) 
        
        Master4ClassificationwithClassAppQuantAug <- # to implement ramp you can manipulate here
          Master4ClassificationwithClass %>% 
          mutate_(.dots = setNames(paste("QuantileFunction(",
                                         names(ThrehsoldsTraining),
                                         ",",
                                         ThrehsoldsTraining,")",
                                         sep = ""),
                                   names(ThrehsoldsTraining)))
        
        rownames(Master4ClassificationwithClassAppQuantAug) <- rownames(Master4ClassificationwithClass)
      }
      #if(Methods[k]=="Logit"){
      
      thislogit <- glmnet(x= Master4ClassificationwithClassAppQuantAug[trainSamples,] %>% 
                            dplyr::select(-CancerStatus) %>% 
                            as.matrix(),
                          y=Master4ClassificationwithClassAppQuantAug[trainSamples,"CancerStatus"],
                          #lambda = lambda.best,
                          penalty.factor = PenaltyFactor,
                          family="binomial",
                          lower.limits=mylower)
      ProbPredictions[testSamples] <- predict(thislogit,
                                              Master4ClassificationwithClassAppQuantAug[testSamples,]%>%
                                                dplyr::select(-CancerStatus) %>% 
                                                as.matrix(),
                                              type="response",
                                              s=lambda.best) %>% 
        unlist()
      
      for( k in Master4ClassificationwithClassAppQuantAug%>% colnames() %>% setdiff("CancerStatus")){
        z <- Master4ClassificationwithClassAppQuantAug[testSamples,] %>%
          dplyr::select(-CancerStatus) %>%
          as.matrix()
        z[,setdiff(colnames(z),k)] <- 0
        
        ProbPredictionsMinusFeature[testSamples,k]<- predict(thislogit,
                                                             z,
                                                             type="response",
                                                             s=lambda.best) %>% 
          unlist()
        
        z <- Master4ClassificationwithClassAppQuantAug[testSamples,] %>%
          dplyr::select(-CancerStatus) %>%
          as.matrix()
        z[,k] <- 0
        
        ProbPredictionsDropFeature[testSamples,k]<- predict(thislogit,
                                                            z,
                                                            type="response",
                                                            s=lambda.best) %>% 
          unlist()
      }
      #}else if(Methods[k]=="RF"){
      
      
      
      
      #  ThisRF <- randomForest(as.factor(CancerStatus)~.,
      #                         data =Master4ClassificationwithClassAppQuantAug[,c("CancerStatus",RFFeatures)] ,
      #                         subset=trainSamples,
      #                         cutoff = c(0.15,0.85))
      #  ProbPredictions[testSamples] <- predict(ThisRF,Master4ClassificationwithClassAppQuantAug[testSamples,],type="prob")[,"TRUE"]
      #}
    }
    
    ThisThresh99 <- ProbPredictions[which(Master4ClassificationwithClass$CancerStatus==F)] %>%
      quantile(SpecSet)
    Allpredictions <- ProbPredictions>ThisThresh99
    AllThresh99NIter[[j]] <- ThisThresh99
    AllpredictionsNIter[[j]] <- Allpredictions
    AllProbPredictionsNIter[[j]] <- ProbPredictions
    AllProbPredictionsMinusFeatureNIter[[j]] <- as.data.frame(ProbPredictionsMinusFeature) %>%
      rownames_to_column("BV") %>%
      mutate( CancerStatus =
                Master4ClassificationwithClass[rownames(ProbPredictionsMinusFeature),
                                               "CancerStatus"]) 
    AllProbPredictionsDropFeatureNIter[[j]] <- as.data.frame(ProbPredictionsDropFeature) %>%
      rownames_to_column("BV") %>%
      mutate( CancerStatus =
                Master4ClassificationwithClass[rownames(ProbPredictionsMinusFeature),
                                               "CancerStatus"])
    AllSensitivity[[j]] <- Allpredictions[which(Master4ClassificationwithClass$CancerStatus==T)]  %>% 
      mean()
    
    Master4ClassificationwithClassAll[[j]] <- cbind(Master4ClassificationwithClass,
                                                    fold=
                                                      MutationsThisRound[rownames(Master4ClassificationwithClass),"fold"])
    
    #} 
    
    
  }
  
  return(list(Master4ClassificationwithClassAll=Master4ClassificationwithClassAll,
              AllThresh99NIter=AllThresh99NIter,
              AllpredictionsNIter=AllpredictionsNIter,
              AllProbPredictionsNIter=AllProbPredictionsNIter,
              AllProbPredictionsMinusFeatureNIter =AllProbPredictionsMinusFeatureNIter,
              AllProbPredictionsDropFeatureNIter = AllProbPredictionsDropFeatureNIter,
              AllSensitivity=AllSensitivity,
              AllSensitivitycvglmnet=max(AllSensity,na.rm = T), # This comes out of Logit cv and more reliable
              AllcvglmnetCalls = Alllogit$fit.preval[,which.max(AllSensity)],
              ProteinsHigherCancer=ProteinsHigherCancer,
              NoMutations=is.null(Mutations),
              Proteins2Use = Proteins2Use,
              Alllogit=Alllogit,
              AlllogitFeatures=AlllogitFeatures,
              #RFFeatures=RFFeatures,
              lambda.best=lambda.best))
  
}



#Preparing Data

Master_Protein_Dataset_12_9_17_Logist <- read.xlsx2("MEGA Protein Dataset for Logistic Regression, December 9, 2017.xlsx",
                                                    sheetIndex = 1,
                                                    stringsAsFactors = F)


Master_Protein_Dataset_12_9_17_RF <- read.xlsx2("MEGA Protein Dataset for Random Forest, December 9, 2017.xlsx",
                                                sheetIndex = 1,
                                                stringsAsFactors = F)

Master_Protein_Filtered_RF <-
  #Master_Protein_Dataset_20_10_17 %>%
  Master_Protein_Dataset_12_9_17_RF %>%
  #filter(!(BV. %in% PostSample)) %>%
  #Master_Protein_Dataset_9_16_17 %>%
  #filter(BV. %in% Master_Protein_Dataset_1_10_17$BV.)%>%
  dplyr::select(-`LRG.1`,-Vitronectin)


TissueType_RFData <- setNames(Master_Protein_Filtered_RF$Sample.Type,Master_Protein_Filtered_RF$BV.)
# ajdusting for limits of detection for each protein?
OnlyProteins_RF <- Master_Protein_Filtered_RF %>%

  dplyr::select(-BV.,-Stage,-Sample.Type) %>% 
  apply(MARGIN = 2,FUN = as.numeric) %>% 
  as.data.frame()


OnlyProteins_RF[is.na(OnlyProteins_RF)] <- 0

colnames(OnlyProteins_RF) <- gsub(colnames(OnlyProteins_RF),pattern = "[.]",replacement = "")
rownames(OnlyProteins_RF) <- Master_Protein_Filtered_RF$BV.


Master_Protein_Filtered_LR <-
  #Master_Protein_Dataset_20_10_17 %>%
  Master_Protein_Dataset_12_9_17_Logist %>%
  #filter(!(BV. %in% PostSample)) %>%
  #Master_Protein_Dataset_9_16_17 %>%
  #filter(BV. %in% Master_Protein_Dataset_1_10_17$BV.)%>%
  dplyr::select(-`LRG.1`,-Vitronectin)

# ajdusting for limits of detection for each protein?
OnlyProteins_LR <- Master_Protein_Filtered_LR %>%
  dplyr::select(-BV.,-Stage,-Sample.Type) %>% 
  apply(MARGIN = 2,FUN = as.numeric) %>% 
  as.data.frame()

OnlyProteins_LR[is.na(OnlyProteins_LR)] <- 0

colnames(OnlyProteins_LR) <- gsub(colnames(OnlyProteins_LR),pattern = "[.]",replacement = "")
rownames(OnlyProteins_LR) <- Master_Protein_Filtered_LR$BV.



## CtdNA

load("maxValuesPerSample_20171209_FORJosh.rda")
sampleTable_ALLModified <- 
  sampleTable_ALL[,c("iteration","fold","maxOmega","CosmicCount","Sample.Category")] %>% 
  rownames_to_column("BV") %>% 
  mutate(BV=ifelse(iteration==1,BV,substr(BV,1,nchar(BV)-1))) 

NewCtDNA_ALLRounds <- sampleTable_ALLModified %>% 
  filter(BV %in% Master_Protein_Filtered_LR$BV.)%>% 
  transmute(BV=BV,
            iteration=iteration,
            fold=fold,
            maxOmega= maxOmega#,maxOmegaThresh),
            #CosmicCount=CosmicCount,
            #maxOmegaCount= QFunclog(CosmicCount*maxOmega,maxOmegaCountThresh)
  ) 

NoNewCtDNA <- setdiff(Master_Protein_Filtered_LR$BV.,NewCtDNA_ALLRounds$BV)
NoNewCtDNA_ALLRounds <- bind_rows(sapply(unique(NewCtDNA_ALLRounds$iteration),
                                         FUN = function(x) 
                                           data.frame(BV=NoNewCtDNA,
                                                      iteration=x,
                                                      fold=as.numeric(balanced.cv.fold(NoNewCtDNA,
                                                                                       length(unique(NewCtDNA_ALLRounds$fold)))),
                                                      maxOmega=0),
                                         simplify = F) )

NewCtDNA_ALLRounds_Modified <- rbind(NewCtDNA_ALLRounds,NoNewCtDNA_ALLRounds)



### RF Analysis#########
#### with CtDNA

OnlyTenProteinsNoCosmicRemovedNewTemplatesRFData <-cv.CancervsNormal.QuantileReplace.Mutatation.knownCV(
  CancerStatus = 
    1*(Master_Protein_Filtered_RF$Sample.Type !="Normal"),
  Proteins = 
    OnlyProteins_RF %>%
    dplyr::select(
      CA199,
      CEA,
      CA125,
      #AFP,
      Prolactin,
      HGF,
      OPN	,
      TIMP1,
      #Follistatin,
      #GCSF,
      #HE4,
      #CA153,
      #IL6,
      #Midkine,
      Myeloperoxidase#,
      #CYFRA211,
      #Galectin3#,
      #Thrombospondin2
    )   %>%
    as.matrix() ,
  
  ProteinsHigherCancer = T,
  NewCtDNA_ALLRounds_Modified %>% 
    mutate(maxOmega=QFuncjumpramp(maxOmega,0)
           #,
           #maxOmegaCount=QFunclog(maxOmegaCount,Inf) 
    ),
  FeaturesToQuantile = 
    c("CA199",
      "CEA",
      "CA125",
      "Prolactin",
      "HGF",
      "OPN"	,
      "TIMP1",
      #"IL6",
      #"Midkine",
      "Myeloperoxidase"),
  QuantileFunction = function(x,Threshold) 
    ifelse(x>Threshold,x,0),
  QuantileF = 0.95,
  ForcePositiveCoeff = F,
  NoLassoPenaly = T )


### Removing Protein Results
mainPath <- "Results/"
dir.create(mainPath,recursive = T)

sink(paste(mainPath,"ResultSummary.txt",sep = ""))
print("---------RF Data Results with 8 Proteins---------")
printRemoveDropResults(ThisResult = OnlyTenProteinsNoCosmicRemovedNewTemplatesRFData)

ProbsRFData <- ProbsEachRound(OnlyTenProteinsNoCosmicRemovedNewTemplatesRFData,
                              csvFile = paste(mainPath,"RFDatawithCtdna8proteins.csv",sep = "") )

LudaObject <- Generate4LudaObject(ThisResult = OnlyTenProteinsNoCosmicRemovedNewTemplatesRFData,
                    rdaFile = paste(mainPath,"ForLuda8proteins.rda",sep = "") )




#Print scores
### Save Logit Scores
dir.create(paste(mainPath,"OnlyThisProtein/",sep = ""),recursive = T)
writeCSVMultiple(rownames(OnlyTenProteinsNoCosmicRemovedNewTemplatesRFData$AllProbPredictionsMinusFeatureNIter) %>% 
                   sapply(FUN = function(x) OnlyTenProteinsNoCosmicRemovedNewTemplatesRFData$AllProbPredictionsMinusFeatureNIter[[x,"LR"]]%>% 
                            mutate(Threshold=OnlyTenProteinsNoCosmicRemovedNewTemplatesRFData$AllThresh99NIter[[x,"LR"]]),simplify = F),
                 file = paste(mainPath,"OnlyThisProtein/",sep = ""))

### Save Logit Scores
dir.create(paste(mainPath,"RemoveThisPrtoein/",sep = ""),recursive = T)

writeCSVMultiple(rownames(OnlyTenProteinsNoCosmicRemovedNewTemplatesRFData$AllProbPredictionsDropFeatureNIter) %>% 
                   sapply(FUN = function(x) OnlyTenProteinsNoCosmicRemovedNewTemplatesRFData$AllProbPredictionsDropFeatureNIter[[x,"LR"]]%>% 
                            mutate(Threshold=OnlyTenProteinsNoCosmicRemovedNewTemplatesRFData$AllThresh99NIter[[x,"LR"]]),simplify = F),
                 file = paste(mainPath,"RemoveThisPrtoein/",sep = ""))


### without CtdNA


OnlyTenProteinsNoCosmicRemovedNewTemplatesRFDatawithoutCtDNA <-cv.CancervsNormal.QuantileReplace.Mutatation.knownCV(
  CancerStatus = 
    1*(Master_Protein_Filtered_RF$Sample.Type !="Normal"),
  Proteins = 
    OnlyProteins_RF %>%
    dplyr::select(
      CA199,
      CEA,
      CA125,
      #AFP,
      Prolactin,
      HGF,
      OPN	,
      TIMP1,
      #Follistatin,
      #GCSF,
      #HE4,
      #CA153,
      #IL6,
      #Midkine,
      Myeloperoxidase#,
      #CYFRA211,
      #Galectin3#,
      #Thrombospondin2
    )   %>%
    as.matrix() ,
  
  ProteinsHigherCancer = T,
  NewCtDNA_ALLRounds_Modified %>% 
    mutate(maxOmega=QFuncjumpramp(maxOmega,Inf)
           #,
           #maxOmegaCount=QFunclog(maxOmegaCount,Inf) 
    ),
  FeaturesToQuantile = 
    c("CA199",
      "CEA",
      "CA125",
      "Prolactin",
      "HGF",
      "OPN"	,
      "TIMP1",
      #"IL6",
      #"Midkine",
      "Myeloperoxidase"),
  QuantileFunction = function(x,Threshold) 
    ifelse(x>Threshold,x,0),
  QuantileF = 0.95,
  ForcePositiveCoeff = F,
  NoLassoPenaly = T )


printRemoveDropResults(ThisResult = OnlyTenProteinsNoCosmicRemovedNewTemplatesRFDatawithoutCtDNA)

ProbsRFData <- ProbsEachRound(OnlyTenProteinsNoCosmicRemovedNewTemplatesRFDatawithoutCtDNA,
                              csvFile = paste(mainPath,"RFDatawithCtdna8proteinswithoutCtdna.csv",sep = "") )




### Removing Each Protein at a time

EightProteins <- c("CA199",
                   "CEA",
                   "CA125",
                   "Prolactin",
                   "HGF",
                   "OPN"	,
                   "TIMP1",
                   "Myeloperoxidase")

SensOnceEightProteins <- setNames(rep(0,length(EightProteins)),EightProteins)

for(i in seq_along(EightProteins)){
  
  EightProteinsRemoved <- setdiff(EightProteins,EightProteins[i])
  
  OnlyTenProteinsNoCosmicRemovedNewTemplatesRFDataRemoveOneProtein <-
    cv.CancervsNormal.QuantileReplace.Mutatation.knownCV(
    CancerStatus = 
      1*(Master_Protein_Filtered_RF$Sample.Type !="Normal"),
    Proteins = 
      OnlyProteins_RF[, EightProteinsRemoved] %>%
      as.matrix() ,
    
    ProteinsHigherCancer = T,
    NewCtDNA_ALLRounds_Modified %>% 
      mutate(maxOmega=QFuncjumpramp(maxOmega,0)
             #,
             #maxOmegaCount=QFunclog(maxOmegaCount,Inf) 
      ),
    FeaturesToQuantile = 
      EightProteinsRemoved,
    QuantileFunction = function(x,Threshold) 
      ifelse(x>Threshold,x,0),
    QuantileF = 0.95,
    ForcePositiveCoeff = F,
    NoLassoPenaly = T )
  
  ProbsRFData <- ProbsEachRound(OnlyTenProteinsNoCosmicRemovedNewTemplatesRFDataRemoveOneProtein,
                                csvFile = paste(mainPath,"RFDatawithCtdna8proteinswithout",EightProteins[i], ".csv",sep = "") )

    
    SensOnceEightProteins[i] <- OnlyTenProteinsNoCosmicRemovedNewTemplatesRFDataRemoveOneProtein$AllSensitivity %>% 
                                                                                        unlist() %>% 
                                                                                          mean()
      
}

cat("--- Remove Protein------\n")
print(t(t(SensOnceEightProteins)))

sink()
