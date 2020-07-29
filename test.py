import cancer_data

import pandas as pd

chunk_iterator = pd.read_csv("temp.csv", chunksize=150000, index_col=0)

for idx, chunk in enumerate(chunk_iterator):

    chunk = chunk.T
    chunk.index = chunk.index.map(lambda x: x.split(".")[0])

    chunk.to_hdf(f"tcga_mx_{idx}.h5", key=f"tcga_mx_{idx}", mode="w", complevel=9, complib="bzip2")
