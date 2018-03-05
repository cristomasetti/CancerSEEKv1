library(openxlsx)
library(plyr)
library(dplyr)

pt = read.xlsx("MEGA PT Data for Lu%2c August 21%2c 2017.xlsx", sheet = 1)
pt$mutID = paste(pt$Chrom, pt$Position,pt$BaseFrom,pt$BaseTo)
pt$Matched.Plasma = as.character(pt$Matched.Plasma)


cancer = c("CRC","Lung","Breast","Pancreas","Ovarian","Esophagus","Liver","Stomach","Small Intestine","Gastric","Ovary","pancreas")

uidThr=200
## threshold for omega value
thres = 1.9

final_cv = data.frame()
fp_cv = data.frame()

for(m in 1:10){
  allResults = c()
  cv = matrix(0, nrow = 10, ncol = 3)
  colnames(cv) = c("FP","Sens", "Con")
  
  for(c in 1:10){
    
    load(paste0("PvalueRatio1209cv_5perc_bl_",m,"_", c,"_1208data.rda"))

    result$fail = (result$UID1<uidThr)+(result$UID2<uidThr)+(result$UID3<uidThr)+
      (result$UID4<uidThr)+(result$UID5<uidThr)+(result$UID6<uidThr)
    result = result[result$fail<=2,]
    
    ## difference between average MAF in the test and the max MAF in the normal controls
    result$diff = result$aveMAF-result$max
    result$diff_r = (result$aveMAF-result$max)/(result$max+10e-6)
    result$PT.Avg.MAF[is.na(result$PT.Avg.MAF)] = -1
    result$PT.Avg.MAF = as.numeric(result$PT.Avg.MAF)
    
    ratio = matrix(c(result$r1, result$r2,result$r3,result$r4, result$r5,result$r6), nrow = nrow(result), ncol = 6)
    uid = matrix(c(result$UID1,result$UID2, result$UID3,result$UID4,result$UID5,result$UID6), nrow = nrow(result), ncol = 6)
    
    order.ratio = t(apply(ratio, 1, order))
    
    for(i in 1:nrow(ratio)){
      ratio[i,] = ratio[i,][order.ratio[i,]]
      uid[i,] = uid[i,][order.ratio[i,]]
    }
    
    uid[uid<uidThr] = 0
    
    ## eliminate the min and max wells in the testSet
    result$omega = rowSums(log(ratio[,-c(1,6)])*uid[,-c(1,6)])/rowSums(uid[,-c(1,6)])

    result$omega[is.na(result$omega)] = Inf
    result$class = FALSE
    result$class = (result$omega>=thres)
    result$iteration = m
    result$fold = c
    ## normal plasmas in the testSet
    nltemp = setdiff(nlpls,nlt)
    
    ## summarizing results
    summ = ddply(result, .(Template,Sample.Category), summarise, PT = max(PT.Avg.MAF*class), class = max(class))
    
    cv[c,1] = sum(summ$class[summ$Template %in% nltemp])/length(nltemp)
    cv[c,2] = sum(summ$class[summ$Sample.Category %in% cancer])/nrow(summ[summ$Sample.Category %in% cancer,])
    
    check = summ[summ$class==TRUE & (summ$Sample.Category %in% cancer) & (summ$Template %in% pt$Matched.Plasma),]
    cv[c,3] = sum(check$PT>=1)/sum(check$PT>=0)
    
    allResults = rbind(allResults,result)
    
  }
  save(allResults, file = paste0("allResults_fromLuMethod_",m,"_20171209.rda"))

  print(m)
  cv = data.frame(cv)
  final_cv = rbind.data.frame(final_cv,cv)
}



sampleTable_ALL = data.frame()
for(m in 1:10){
  load(file = paste0("allResults_fromLuMethod_",m,"_20171209.rda"))
  
  # subset by sample omega values and mutations
  sampleOmegaList = tapply(allResults[,'omega'], INDEX = allResults[,'Template'], FUN = function(x)return(x))
  sampleMutList = tapply(allResults[,'mutID'], INDEX = allResults[,'Template'], FUN = function(x)return(x))
  sampleCosmicList = tapply(allResults[,'CosmicCount'], INDEX = allResults[,'Template'], FUN = function(x)return(x))
  sampleGeneList = tapply(allResults[,'Gene'], INDEX = allResults[,'Template'], FUN = function(x)return(x))
  sampleAmpList = tapply(allResults[,"ampMatchName"], INDEX = allResults[,'Template'], FUN = function(x)return(x))
  
  # find max value
  sampleMaxS = unlist(tapply(allResults[,'omega'], INDEX = allResults[,'Template'], FUN = which.max, simplify = T))
  
  maxS = sapply(names(sampleMaxS), function(i) sampleOmegaList[[i]][sampleMaxS[i]])
  maxSmut = sapply(names(sampleMaxS), function(i) sampleMutList[[i]][sampleMaxS[i]])
  cosm = sapply(names(sampleMaxS), function(i) sampleCosmicList[[i]][sampleMaxS[i]])
  gene = sapply(names(sampleMaxS), function(i) sampleGeneList[[i]][sampleMaxS[i]])
  amp = sapply(names(sampleMaxS), function(i) sampleAmpList[[i]][sampleMaxS[i]])
  
  
  sampleTable = allResults[,c('Template','Sample.Category',"iteration","fold")]
  sampleTable = sampleTable[!duplicated(sampleTable[,1]),]
  rownames(sampleTable) = sampleTable[,1]
  sampleTable = cbind(sampleTable,maxOmega = maxS[rownames(sampleTable)], maxO_mut = maxSmut[rownames(sampleTable)],
                      CosmicCount = cosm[rownames(sampleTable)], gene = gene[rownames(sampleTable)], ampMatchName = amp[rownames(sampleTable)])
  
  sampleTable$remove = paste(sampleTable$Template,sampleTable$maxO_mut)
  sampleTable = sampleTable[,-1]
  
  uid_matrix = allResults[,c(1:51)]
  uid_matrix_subset = uid_matrix[match(sampleTable$remove,uid_matrix$remove),]
  
  sampleTable = cbind(sampleTable,uid_matrix_subset)
  
  #save(sampleTable, file = 'omegaPerSample_20171019.rda')
  
  # add maximum diff and diff_r per sample
  diff = unlist(tapply(allResults[,'diff'], INDEX = allResults[,'Template'], FUN = max, na.rm = T, simplify = T))
  diff_r = unlist(tapply(allResults[,'diff_r'], INDEX = allResults[,'Template'], FUN = max, na.rm = T, simplify = T))
  sampleTable_diff = cbind(sampleTable,maxDiff = diff[rownames(sampleTable)],maxDiff_r = diff_r[rownames(sampleTable)])
  sampleTable_ALL = rbind.data.frame(sampleTable_ALL, sampleTable_diff)
}
save(sampleTable_ALL, file = 'maxValuesPerSample_20171209_FORJosh.rda')


