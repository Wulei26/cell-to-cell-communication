library(Seurat)
library(tidyverse)
library(data.table)
library(readr)
library(qs)
crc = qread("/public2022/wulei/research/counts/tumor/CRC.qs")

crc = subset(crc, dataset == 8)
crc = CreateSeuratObject(crc[["RNA"]]@counts,
                        assay = "RNA", min.cells = 3,
                        meta.data = crc@meta.data)
## 保留  Epithelial Fibroblast Myeloid
crc = subset(crc, (major %in% c("Epithelial", "Fibroblast", "Myeloid")) | (revised_subtype == "Treg"))

crc$cell_type_need = ifelse(crc$revised_subtype %in% c("Cycling TAM","FOLR2+TAM","FTL+TAM","SPP1+APOE+TAM"),"TAM",
                        ifelse(crc$revised_subtype %in% c("AP CAF","detox-iCAF","ecm-myCAF","Epi-like CAF","IFNg-iCAF","IL-iCAF","wound-myCAF"),"CAF",
                        ifelse(crc$revised_subtype == "Treg","Treg",
                        ifelse(crc$major == "Epithelial",crc$revised_subtype,"NA"))))

crc = subset(crc, cell_type_need != "NA")
set.seed(10086)
## 手动进行一个 downsample
Idents(crc) = "cell_type_need"
crcsub = subset(crc, downsample = 5000)
crcsub = NormalizeData(crcsub)

library(SeuratDisk)
# library(Seurat)      
#seurat2h5seurat中间过渡	
SaveH5Seurat(crcsub,filename="/public2022/wulei/research/CCC/CRC/raw_data/merge/crcsub.h5seurat", overwrite = TRUE)
#数据转为最终h5ad格式
Convert("/public2022/wulei/research/CCC/CRC/raw_data/merge/crcsub.h5seurat", dest = "/public2022/wulei/research/CCC/CRC/raw_data/merge/crcsub.h5ad", overwrite = TRUE)


### 分样本写入
file_path = "/public2022/wulei/research/CCC/CRC/raw_data/data/"
for (pi in as.character(unique(crc$patientid))) {
   cat(pi)
   scesub = subset(crc, patientid == pi)
   scesub = NormalizeData(scesub)
   ## 存储data
   ##先创建一个文件夹
   dir.create(paste0(file_path,pi))
   data = as.data.frame(scesub@assays$RNA@data) %>% rownames_to_column(var = "Gene")
   row_sums_greater_than_zero <- apply(data[, 2:ncol(data)], 1, function(x) sum(x) > 0)

   # 使用逻辑向量筛选数据框
   data <- data[row_sums_greater_than_zero, ]

   fwrite(data, paste0(file_path,pi,"/",pi,"_Normalized_data.txt"), row.names=F, sep='\t')

   ## 写入metadata
   meta <- data.frame(Cell=rownames(scesub@meta.data),
                       cell_type=scesub@meta.data$cell_type_need)
     
   fwrite(meta, paste0(file_path,pi,"/",pi,"_meta.txt"), row.names=F, sep='\t')

   qsave(scesub,file = paste0(file_path,pi,"/",pi,".qs"))
}
