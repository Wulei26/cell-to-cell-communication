import pandas as pd
import sys
import os
print(sys.version)
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method


raw_data = '/public2022/wulei/research/CCC/CRC/raw_data/data'
base_dir = [os.path.join(raw_data, file) for file in os.listdir(raw_data)]
base_dir.remove('/public2022/wulei/research/CCC/CRC/raw_data/data/preparedata.R')

cpdb_file_path = '/public2022/wulei/research/CCC/cpdb_tutorial/v5.0.0/cellphonedb.zip'


for dir in base_dir:
    ## 创建一个结果文件夹
    os.mkdir(dir + "/result")
    sample_name = dir[-4: ]
    counts_file_path = dir + "/" + sample_name + "_Normalized_data.txt"
    meta_file_path = dir + "/" + sample_name + "_meta.txt"
    out_path = dir + "/result"
    if os.path.exists(out_path):
        try:
            cpdb_results = cpdb_statistical_analysis_method.call(
            cpdb_file_path = cpdb_file_path,                 # mandatory: CellphoneDB database zip file.
            meta_file_path = meta_file_path,                 # mandatory: tsv file defining barcodes to cell label.
            counts_file_path = counts_file_path,             # mandatory: normalized count matrix - a path to the counts file, or an in-memory AnnData object
            counts_data = 'hgnc_symbol',                     # defines the gene annotation in counts matrix.
            score_interactions = True,                       # optional: whether to score interactions or not. 
            output_path =  out_path,                          # Path to save results.
            output_suffix = sample_name                             # Replaces the timestamp in the output files by a user defined string in the  (default: None).
            )
        except OSError as e:
            print(f"目录不存在: {e}")
    else:
        print(f"目录 {out_path} 不存在")   
