library(openxlsx)
library(dplyr)
## prepare data
directory_path = "/public2022/wulei/research/CCC/CRC/raw_data/data"
sampleNames = folders <- list.dirs(directory_path, full.names = FALSE, recursive = FALSE)

setwd("/public2022/wulei/research/CCC/CRC/raw_data/data")

tumor_subtypes = c("BEST2+Goblet","BEST4+Enterocyte","Colonocyte","Cycling Epi","Enterocyte",
                "Enteroendocrine","LGR5+ stem cell","Microfold","Paneth","Stem cells","Tuft cell")
anames <- c("id_cp_interaction", "interacting_pair", "partner_a", "partner_b",
    "gene_a", "gene_b", "secreted", "receptor_a", "receptor_b", "annotation_strategy",
    "is_integrin")
## Treg
tregToTumorCell = paste0("Treg|",tumor_subtypes)
tumorCellToTreg = paste0(tumor_subtypes,"|Treg")

getUniqueAllId = function(dat) {
    idMap = c()
    for(i in sampleNames) {
        vec = dat[[i]]
        vec = vec[!is.na(vec)]
        vec = unique(vec)
        idMap = c(idMap,vec )
    }
    idMap = unique(idMap)
    if("noValue" %in% idMap) {
        idMap = idMap[idMap != "noValue"]
    }
    return(idMap)
}

resFinal = list()
for (sub in tregToTumorCell) {
    res = data.frame(matrix(data = NA,ncol = length(sampleNames),nrow = 500))
    colnames(res) = sampleNames
    for (sam in sampleNames) {
        sig = read.table(paste0("./",sam,"/result/","statistical_analysis_significant_means_",sam,".txt"),sep = "\t",check.names = FALSE,header = TRUE)
        if (sub %in% colnames(sig)) {
            dat.sam = sig[!is.na(sig[[sub]]),]
            if(dim(dat.sam)[1] == 0) {
                res[[sam]] = c("noValue",rep(NA,499))
            } else {
                res[[sam]] = c(dat.sam$id_cp_interaction,rep(NA,500 - length(dat.sam$id_cp_interaction)))
            }
        } else {
            res[[sam]] = c("noValue",rep(NA,499))
        }
    }
    resFinal[[sub]] = res
}

for(i in names(resFinal)) {
    DF = resFinal[[i]]
    cat(i,dim(DF)[2],"\n")
    # cat(cols_all_na,"\n")
}

## 获取 "Treg|BEST2+Goblet"的 所有存在的通讯id
resSearch = list()
for(sm in names(resFinal)) {
    df = resFinal[[sm]]
    nid = getUniqueAllId(df)
    datMaster = data.frame(matrix(data = "NA",ncol = length(sampleNames), nrow = length(nid)))
    colnames(datMaster) = sampleNames
    rownames(datMaster) = nid
    for(pid in sampleNames) {
        pidVec = df[[pid]]
        pidVec = unique(pidVec)
        pidVec = pidVec[!is.na(pidVec)]
        for(id in pidVec) {
            if(id %in% nid) {
                datMaster[id,pid] = 1
            }
        }       
    }
    datMaster[datMaster != 1] = 0
    count_ones <- apply(datMaster, 1, function(x) sum(x == 1))
    datMaster$CountOfOnes <- count_ones
    datMaster = datMaster[order(datMaster$CountOfOnes,decreasing = TRUE),]
    resSearch[[sm]] = datMaster    
}

resFinal2 = list()
for (sub in tumorCellToTreg) {
    res = data.frame(matrix(data = NA,ncol = length(sampleNames),nrow = 500))
    colnames(res) = sampleNames
    for (sam in sampleNames) {
        sig = read.table(paste0("./",sam,"/result/","statistical_analysis_significant_means_",sam,".txt"),sep = "\t",check.names = FALSE,header = TRUE)
        if (sub %in% colnames(sig)) {
            dat.sam = sig[!is.na(sig[[sub]]),]
            if(dim(dat.sam)[1] == 0) {
                res[[sam]] = c("noValue",rep(NA,499))
            } else {
                res[[sam]] = c(dat.sam$id_cp_interaction,rep(NA,500 - length(dat.sam$id_cp_interaction)))
            }
        } else {
            res[[sam]] = c("noValue",rep(NA,499))
        }
    }
    resFinal2[[sub]] = res
}
# for(i in names(resFinal2)) {
#     DF = resFinal2[[i]]
#     cat(i,dim(DF)[2],"\n")
#     # cat(cols_all_na,"\n")
# }

## 获取 "Treg|BEST2+Goblet"的 所有存在的通讯id
resSearch2 = list()
for(sm in names(resFinal2)) {
    df = resFinal2[[sm]]
    nid = getUniqueAllId(df)
    datMaster = data.frame(matrix(data = "NA",ncol = length(sampleNames), nrow = length(nid)))
    colnames(datMaster) = sampleNames
    rownames(datMaster) = nid
    for(pid in sampleNames) {
        pidVec = df[[pid]]
        pidVec = unique(pidVec)
        pidVec = pidVec[!is.na(pidVec)]
        for(id in pidVec) {
            if(id %in% nid) {
                datMaster[id,pid] = 1
            }
        }       
    }
    datMaster[datMaster != 1] = 0
    count_ones <- apply(datMaster, 1, function(x) sum(x == 1))
    datMaster$CountOfOnes <- count_ones
    datMaster = datMaster[order(datMaster$CountOfOnes,decreasing = TRUE),]
    resSearch2[[sm]] = datMaster    
}

#  对 resSearch 进行处理
resSearchEnd = list()
for(nas in names(resSearch)) {
    searchSub = resSearch[[nas]]
    searchSubCopy = searchSub 
    searchSub$id_cp_interaction = rownames(searchSub)
    searchSub = searchSub[, c(ncol(searchSub), 1:(ncol(searchSub) - 1))]
    rownames(searchSub) = NULL
    for(idInter in searchSub$id_cp_interaction) {
        # idDat = subset(allResult, id_cp_interaction == idInter, select = c(anames,nas,"patientid"))
        for(pid in sampleNames) {
            sig = read.table(paste0("./",pid,"/result/","statistical_analysis_significant_means_",pid,".txt"),sep = "\t", check.names = FALSE, header = TRUE)
            if(searchSub[searchSub$id_cp_interaction == idInter, pid] == 1){
                searchSubCopy[idInter,pid] = subset(sig, id_cp_interaction == idInter, select = nas)
            } else {
                searchSubCopy[idInter,pid] = NA
            }
        } 
    }
    resSearchEnd[[nas]] = searchSubCopy
}

#  对 resSearch2 进行处理
resSearchEnd2 = list()
for(nas in names(resSearch2)) {
    searchSub = resSearch2[[nas]]
    searchSubCopy = searchSub 
    searchSub$id_cp_interaction = rownames(searchSub)
    searchSub = searchSub[, c(ncol(searchSub), 1:(ncol(searchSub) - 1))]
    rownames(searchSub) = NULL
    for(idInter in searchSub$id_cp_interaction) {
        # idDat = subset(allResult, id_cp_interaction == idInter, select = c(anames,nas,"patientid"))
        for(pid in sampleNames) {
            sig = read.table(paste0("./",pid,"/result/","statistical_analysis_significant_means_",pid,".txt"),sep = "\t", check.names = FALSE, header = TRUE)
            if(searchSub[searchSub$id_cp_interaction == idInter, pid] == 1){
                searchSubCopy[idInter,pid] = subset(sig, id_cp_interaction == idInter, select = nas)
            } else {
                searchSubCopy[idInter,pid] = NA
            }
        } 
    }
    resSearchEnd2[[nas]] = searchSubCopy
}


## 将某个id 的 anames 信息提取出来

#  将所有的结果进行合并
for(sam in sampleNames) {
    sig = read.table(paste0("./",sam,"/result/","statistical_analysis_significant_means_",sam,".txt"),sep = "\t",check.names = FALSE,header = TRUE)
    sig = subset(sig, select = anames)
    sig$patientid = sam
    if(sam == sampleNames[1]) {
        allResult = sig
    }else {
        allResult = rbind(allResult,sig)
    }
}

## 对resSearchEnd 进行处理
resSearchEndResult = list()
for(sub in names(resSearchEnd)) {
    resSearchEndSub = resSearchEnd[[sub]]
    ## 将其行名转化为id_cp_interaction这一列
    resSearchEndSub$id_cp_interaction = rownames(resSearchEndSub)
    # rownames(resSearchEndSub) == NULL
    allResultsub = subset(allResult, allResult$id_cp_interaction %in% rownames(resSearchEndSub))
    allResultsub = subset(allResultsub, select = -patientid)
    rownames(allResultsub) = NULL
    allResultsub = allResultsub[!duplicated(allResultsub),]
    if(dim(resSearchEndSub)[1] == dim(allResultsub)[1]) {
        ##  可以进行合并
        mergeRes = right_join(allResultsub, resSearchEndSub, by = "id_cp_interaction")
        mergeRes = mergeRes[order(mergeRes$CountOfOnes, decreasing = TRUE),]
        rownames(mergeRes) = NULL
    } else {
        cat(sub, "\t Somthing Wrong ！！！","\n")
    }
    resSearchEndResult[[sub]] = mergeRes
}

## 对resSearchEnd2 进行处理
resSearchEndResult2 = list()
for(sub in names(resSearchEnd2)) {
    resSearchEndSub = resSearchEnd2[[sub]]
    ## 将其行名转化为id_cp_interaction这一列
    resSearchEndSub$id_cp_interaction = rownames(resSearchEndSub)
    # rownames(resSearchEndSub) == NULL
    allResultsub = subset(allResult, allResult$id_cp_interaction %in% rownames(resSearchEndSub))
    allResultsub = subset(allResultsub, select = -patientid)
    rownames(allResultsub) = NULL
    allResultsub = allResultsub[!duplicated(allResultsub),]
    if(dim(resSearchEndSub)[1] == dim(allResultsub)[1]) {
        ##  可以进行合并
        mergeRes = right_join(allResultsub, resSearchEndSub, by = "id_cp_interaction")
        mergeRes = mergeRes[order(mergeRes$CountOfOnes, decreasing = TRUE),]
        rownames(mergeRes) = NULL
    } else {
        cat(sub, "\t Somthing Wrong ！！！","\n")
    }
    resSearchEndResult2[[sub]] = mergeRes
}
# resSearchEndResult 和 resSearchEndResult2 进行保存

write.xlsx(x = resSearchEndResult, file = "/public2022/wulei/research/CCC/CRC/raw_data/merger_split_result/Treg/TregToTumorSubtypes.xlsx", rowNames = FALSE, colNames = TRUE, names = names(resSearchEndResult))

write.xlsx(x = resSearchEndResult2, file = "/public2022/wulei/research/CCC/CRC/raw_data/merger_split_result/Treg/TumorSubtypesToTreg.xlsx", rowNames = FALSE, colNames = TRUE, names = names(resSearchEndResult2))
