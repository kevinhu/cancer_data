---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.2'
      jupytext_version: 1.4.2
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
sys.path.append(os.path.relpath("../helper"))

import name_mappings

```

# CRISPR knockdown scores

```python
sanger_model_info = pd.read_csv("raw/sanger/model_list_latest.csv")
sanger_model_map = dict(zip(sanger_model_info["model_name"],sanger_model_info["BROAD_ID"]))

```

```python
score = pd.read_csv("raw/sanger/EssentialityMatrices/01_corrected_logFCs.tsv",sep="\t",index_col=0).T
score.head(5)
score = score.drop("HT-29v1.0",axis=0)
score = score.drop("HT-29v1.1",axis=0)

score.index = score.index.map(lambda x: sanger_model_map[x])
score.dtype = np.float16

score.to_hdf("processed/sanger/01_corrected_logFCs.h5", key="score")

```

```python
score_bayes = pd.read_csv("raw/sanger/EssentialityMatrices/03_scaledBayesianFactors.tsv",sep="\t",index_col=0).T
score_bayes.head(5)
score_bayes = score_bayes.drop("HT-29v1.0",axis=0)
score_bayes = score_bayes.drop("HT-29v1.1",axis=0)

score_bayes.index = score_bayes.index.map(lambda x: sanger_model_map[x])
score_bayes.dtype = np.float16

score_bayes.to_hdf("processed/sanger/03_scaledBayesianFactors.h5", key="score_bayes")

```
