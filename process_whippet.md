---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.2'
      jupytext_version: 1.5.0
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

```python
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

from collections import Counter, defaultdict

import sys
import os
import re
from os.path import isfile, join

sys.path.append(os.path.relpath("../helper"))


import name_mappings

```

```python
sra_info = pd.read_csv("../data/raw/ccle-whippet/ccle_sra.txt",sep="\t")
run_name_map = dict(zip(sra_info["run_accession"],sra_info["sample_title"]))

```

```python
whippet_path = "../data/raw/ccle-whippet/whippet/"

for sra_idx, sra_id in enumerate(sra_info["run_accession"]):

    if sra_idx % 10 == 0:
        print(sra_idx, end=" ")

    #     if not isfile(whippet_path+sra_id+".gene.tpm.gz"):
    #         print(sra_id)

    ccle_name = run_name_map[sra_id]

    if isfile(whippet_path+sra_id+".gene.tpm.gz"):
        psis = pd.read_csv(whippet_path+sra_id+".psi.gz",
                           compression='gzip', sep='\t')

        psis = psis.dropna(subset=["Psi"])
        psis["id"] = psis["Gene"] + "_" + \
            psis["Type"].astype(str) + "_" + psis["Coord"]
        psis = psis.set_index("id")

        psis = psis[["Psi", "CI_Width", "Complexity", "Entropy", "Strand"]]

        psis["Psi"] = psis["Psi"].astype(np.float16)
        psis["CI_Width"] = psis["CI_Width"].astype(np.float16)
        psis["Complexity"] = psis["Complexity"].astype(str)
        psis["Entropy"] = psis["Entropy"].astype(np.float16)
        psis["Strand"] = psis["Strand"].astype(str)

        psis.to_hdf("../data/intermediate/ccle-whippet/psis/" +
                    sra_id+".hdf", key="psis", complevel=4)

        gene_tpms = pd.read_csv(whippet_path+sra_id+".gene.tpm.gz",
                                compression='gzip', sep='\t', skiprows=1,
                                names=["gene_id", "tpm", "counts", "x"]).iloc[:, :3]
        
        gene_tpms = gene_tpms.set_index("gene_id")

        gene_tpms.to_hdf("../data/intermediate/ccle-whippet/gene_tpms/" +
                         sra_id+".hdf", key="gene_tpms", complevel=4)

        isoform_tpms = pd.read_csv(
            whippet_path+sra_id+".isoform.tpm.gz",
            compression='gzip', sep='\t', skiprows=1,
            names=["transcript_id", "tpm", "counts", "x"]).iloc[:, :3]
        
        isoform_tpms = isoform_tpms.set_index("transcript_id")

        isoform_tpms.to_hdf("../data/intermediate/ccle-whippet/isoform_tpms/" +
                            sra_id+".hdf", key="isoform_tpms", complevel=4)
```

```python
merged_psis = []
merged_gene_tpms = []
merged_isoform_tpms = []

for sra_idx, sra_id in enumerate(sra_info["run_accession"]):

    if sra_idx % 10 == 0:
        print(sra_idx, end=" ")

    ccle_name = run_name_map[sra_id]

    if isfile(whippet_path+sra_id+".gene.tpm.gz"):
        
        ccle_name = run_name_map[sra_id]
        
        psis = pd.read_hdf("../data/intermediate/ccle-whippet/psis/" +
                    sra_id+".hdf", key="psis")
        
        psis = psis[psis["CI_Width"] < 0.33]
        
        merged_psis.append(pd.Series(psis["Psi"],name=ccle_name))
        
        gene_tpms = pd.read_hdf("../data/intermediate/ccle-whippet/gene_tpms/" +
                         sra_id+".hdf", key="gene_tpms")
        merged_gene_tpms.append(pd.Series(gene_tpms["tpm"],name=ccle_name))
        
        isoform_tpms = pd.read_hdf("../data/intermediate/ccle-whippet/isoform_tpms/" +
                            sra_id+".hdf", key="isoform_tpms")
        merged_isoform_tpms.append(pd.Series(isoform_tpms["tpm"],name=ccle_name))
```

```python
merged_psis = pd.concat(merged_psis,axis=1,join="outer",sort=True).T
merged_gene_tpms = pd.concat(merged_gene_tpms,axis=1,join="outer",sort=True).T
merged_isoform_tpms = pd.concat(merged_isoform_tpms,axis=1,join="outer",sort=True).T

merged_psis.index = merged_psis.index.map(lambda x: name_mappings.name_map[x])
merged_gene_tpms.index = merged_gene_tpms.index.map(lambda x: name_mappings.name_map[x])
merged_isoform_tpms.index = merged_isoform_tpms.index.map(lambda x: name_mappings.name_map[x])
```

```python
g19_definitions = pd.read_csv("../data/raw/ccle-whippet/gencode_hg19.v25.tsl1.gtf",sep="\t",skiprows=6,
                              names=["chrom","source","type","start","end",".1","strand",".2","info"])

g19_definitions["ensembl_gene_id"] = g19_definitions["info"].apply(
    lambda x: x.split(";")[0][9:-1])
g19_definitions["gene_name"] = g19_definitions["info"].apply(
    lambda x: x.split("gene_name")[1].split(";")[0][2:-1])
g19_definitions["ensembl_tx_id"] = g19_definitions["info"].apply(
    lambda x: x.split("transcript_id")[1].split(";")[0][2:-1])

# g19_gene_definitions = g19_definitions[g19_definitions["type"]=="gene"]

ensembl_id_map = dict(zip(g19_definitions["ensembl_gene_id"],g19_definitions["gene_name"]))
ensembl_id_map = defaultdict(str, ensembl_id_map)

ensembl_transcript_map = dict(zip(g19_definitions["ensembl_tx_id"],g19_definitions["gene_name"]))
ensembl_transcript_map = defaultdict(str, ensembl_transcript_map)

```

```python
merged_gene_tpms = np.log2(merged_gene_tpms+1)
merged_isoform_tpms = np.log2(merged_isoform_tpms+1)

merged_psis.columns = [ensembl_id_map["_".join(x.split("_")[:2])]+"_"+x for x in merged_psis.columns]
merged_gene_tpms.columns = [ensembl_id_map[x]+"_"+x for x in merged_gene_tpms.columns]
merged_isoform_tpms.columns = [ensembl_transcript_map[x]+"_"+x for x in merged_isoform_tpms.columns]
```

```python
psi_stdevs = merged_psis.std(axis=0)
sns.distplot(psi_stdevs.dropna())
merged_psis = merged_psis[merged_psis.columns[psi_stdevs>0.05]]
```

```python
merged_psis.to_hdf("../data/intermediate/ccle-whippet/merged_psis.hdf",key="merged_psis")
merged_gene_tpms.to_hdf("../data/intermediate/ccle-whippet/merged_gene_tpms.hdf",key="merged_gene_tpms")
merged_isoform_tpms.to_hdf("../data/intermediate/ccle-whippet/merged_isoform_tpms.hdf",key="merged_isoform_tpms")
```
