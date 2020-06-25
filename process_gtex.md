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
from functools import reduce

from scipy.stats import zscore

from collections import Counter, defaultdict

import sys
import os
import re
sys.path.append(os.path.relpath("../helper"))

sys.path.append(os.path.relpath("../huygens"))
sys.path.append(os.path.relpath("../galileo"))

import galileo as gal
import huygens as huy

import name_mappings

```

```python
def concat_cols(df, cols, delim):
    cols_str = [df[x].astype(str) for x in cols]

    return reduce(lambda a, b: a + delim + b, cols_str)
```

# RNAseq


## Manifest

```python
gtex_manifest_1 = pd.read_csv("raw/gtex/E-MTAB-5214.sdrf.txt", sep="\t")
gtex_manifest_2 = pd.read_csv("raw/gtex/E-MTAB-2919.sdrf.txt", sep="\t")

gtex_manifest = pd.concat([
    gtex_manifest_1,
    gtex_manifest_2
], 
    axis=0, sort=True)

sra_gtex_map = dict(
    zip(gtex_manifest["Comment[ENA_RUN]"], gtex_manifest["Source Name"]))
```

## Gene expression

```python
gtex_genex = pd.read_csv("raw/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz",skiprows=2,index_col=0,sep="\t")
```

```python
gtex_genex.index = gtex_genex["Description"] + "_" + gtex_genex.index
gtex_genex.drop(["Description"],axis=1,inplace=True)
gtex_genex = gtex_genex.T
gtex_genex = np.log10(gtex_genex+1)
gtex_genex = gtex_genex.astype(np.float32)
```

```python
gtex_genex.to_hdf("processed/gtex/gtex_genex.hdf",key="gtex_genex",mode="w")
```

## 3' splicing

```python
gtex_a3ss = pd.read_csv("raw/gtex/gtex_merge_graphs_alt_3prime_C2.confirmed.txt.gz",sep="\t")
```

```python
gtex_a3ss["exon_id"] = concat_cols(gtex_a3ss, ["gene_name",
                                           "event_type",
                                           "event_chr",
                                           "event_coordinates",
                                           "alt_region_coordinates"],
                                   "_"
                                  )
    
gtex_a3ss = gtex_a3ss.drop(["event_id",
                    "event_type",
                    "event_chr",
                    "event_coordinates",
                    "alt_region_coordinates",
                    "gene_name"
                   ],axis=1)

gtex_a3ss = gtex_a3ss.set_index("exon_id")
gtex_a3ss.columns = [sra_gtex_map[x[:-8]] for x in gtex_a3ss.columns]
gtex_a3ss = gtex_a3ss.T
gtex_a3ss = gtex_a3ss.astype(np.float32)
```

```python
gtex_a3ss.to_hdf("processed/gtex/gtex_a3ss.hdf",key="gtex_a3ss",mode="w")
```

## Exon skip splicing

```python
gtex_se = pd.read_csv("raw/gtex/gtex_merge_graphs_exon_skip_C2.confirmed.txt.gz",sep="\t")
```

```python
gtex_se["exon_id"] = concat_cols(gtex_se, ["gene_name",
                                           "event_type",
                                           "event_chr",
                                           "event_coordinates",
                                           "alt_region_coordinates"],
                                   "_"
                                  )
    
gtex_se = gtex_se.drop(["event_id",
                    "event_type",
                    "event_chr",
                    "event_coordinates",
                    "alt_region_coordinates",
                    "gene_name"
                   ],axis=1)

gtex_se = gtex_se.set_index("exon_id")
gtex_se.columns = [sra_gtex_map[x[:-8]] for x in gtex_se.columns]
gtex_se = gtex_se.T
gtex_se = gtex_se.astype(np.float32)
```

```python
gtex_se.to_hdf("processed/gtex/gtex_se.hdf",key="gtex_se",mode="w")
```
