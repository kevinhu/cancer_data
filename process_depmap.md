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
import pandas as pd

from collections import Counter, defaultdict

import sys
import os
import re
```

```python
def format_name(x):
    x_split = x.split(" (")
    return x_split[0] + "_" + x_split[1][:-1]
```

# Avana

```python
avana = pd.read_csv("raw/depmap/Achilles_gene_effect_avana_19q4.csv",index_col=0).T

avana.index = avana.index.map(lambda x: format_name(x))
avana = avana.T

avana.to_hdf("processed/depmap/avana.hdf",key="avana")
```

# DRIVE

```python
drive = pd.read_csv("raw/depmap/demeter2-drive_v12-gene-effect.csv",index_col=0).T

drive.index = drive.index.map(lambda x: format_name(x))
drive = drive.T
drive = drive.drop(["GISTT1_GASTROINTESTINAL_TRACT"],axis=0)
drive.index = drive.index.map(lambda x: name_mappings.name_map[x])

drive.to_hdf("processed/depmap/demeter2-drive_v12-gene-effect.hdf",key="drive")
```

# Achilles

```python
achilles = pd.read_csv("raw/depmap/demeter2-achilles_v13-gene-effect.csv",index_col=0).T

achilles.index = achilles.index.map(lambda x: format_name(x))
achilles = achilles.T
achilles.index = achilles.index.map(lambda x: name_mappings.name_map[x])

achilles.to_hdf("processed/depmap/demeter2-achilles_v13-gene-effect.hdf",key="achilles")
```

# Mutation classes


## Damaging

```python
damaging_muts = pd.read_csv("raw/depmap/depmap-mutation-calls_v11-damaging-mutation.csv",index_col=0).T

damaging_muts.index = damaging_muts.index.map(lambda x: format_name(x))
damaging_muts = damaging_muts.T

damaging_muts.to_hdf("processed/depmap/depmap-mutation-calls_v11-damaging-mutation.hdf",key="damaging_muts")
```

## Hotspot

```python
hs_muts = pd.read_csv("raw/depmap/depmap-mutation-calls_v11-hotspot-mutation.csv",index_col=0).T

hs_muts.index = hs_muts.index.map(lambda x: format_name(x))
hs_muts = hs_muts.T

hs_muts.to_hdf("processed/depmap/depmap-mutation-calls_v11-hotspot-mutation.hdf",key="hs_muts")
```

## Other

```python
other_muts = pd.read_csv("raw/depmap/depmap-mutation-calls_v11-other-mutation.csv",index_col=0).T

other_muts.index = other_muts.index.map(lambda x: format_name(x))
other_muts = other_muts.T

other_muts.to_hdf("processed/depmap/depmap-mutation-calls_v11-other-mutation.hdf",key="other_muts")
```

# Drug responses (GDSC)


## Dose-responses

```python
sanger_response = pd.read_csv("raw/depmap/sanger-dose-response.csv")

sanger_response["drug_id"] = sanger_response["DRUG_NAME"].fillna("UNNAMED") + "_" + sanger_response["DRUG_ID"].astype(str) + "_" + sanger_response["DATASET"]

```

```python
sanger_ic50 = pd.pivot_table(sanger_response, 
                             values="IC50_PUBLISHED", 
                             index=["ARXSPAN_ID"],
                             columns="drug_id", 
                             fill_value=np.nan)

sanger_auc = pd.pivot_table(sanger_response, 
                            values="AUC_PUBLISHED",
                            index=["ARXSPAN_ID"], 
                            columns="drug_id", 
                            fill_value=np.nan)

sanger_ic50.to_hdf("processed/depmap/sanger_ic50.h5",key="sanger_ic50",mode="w")
sanger_auc.to_hdf("processed/depmap/sanger_auc.h5",key="sanger_auc",mode="w")
```

# Drug responses (PRISM)


## Logfolds

```python
primary_logfold = pd.read_csv("raw/depmap/primary-screen-public-tentative_v10-primary-replicate-collapsed-logfold-change.csv",index_col=0)
secondary_logfold = pd.read_csv("raw/depmap/secondary-screen-replicate-collapsed-logfold-change.csv",index_col=0)

primary_info = pd.read_csv("raw/depmap/primary-screen-public-tentative_v10-primary-replicate-collapsed-treatment-info.csv")
secondary_info = pd.read_csv("raw/depmap/secondary-screen-replicate-collapsed-treatment-info.csv")

primary_info["format_name"] = primary_info["name"].fillna("UNNAMED") + "_" + primary_info["column_name"]
secondary_info["format_name"] = secondary_info["name"].fillna("UNNAMED") + "_" + secondary_info["column_name"]

primary_name_map = dict(zip(primary_info["column_name"],primary_info["format_name"]))
secondary_name_map = dict(zip(secondary_info["column_name"],secondary_info["format_name"]))

primary_logfold.columns = [primary_name_map[x] for x in primary_logfold.columns]
secondary_logfold.columns = [secondary_name_map[x] for x in secondary_logfold.columns]

primary_logfold = primary_logfold.astype(np.float32)
secondary_logfold = secondary_logfold.astype(np.float32)
```

```python
primary_logfold.to_hdf(
    "processed/depmap/primary_logfold.h5",
    key="primary_logfold",
    mode="w"
)

secondary_logfold.to_hdf(
    "processed/depmap/secondary_logfold.h5",
    key="secondary_logfold",
    mode="w"
)
```

# Copy number

```python
copynumber = pd.read_csv("raw/depmap/CCLE_gene_cn_19q4_public.csv",index_col=0)

copynumber.to_hdf("processed/depmap/CCLE_gene_cn_19q4_public.hdf",key="copynumber")
```

# MSI

```python
is_msi = pd.read_csv("raw/depmap/CCLE_MSI.csv",index_col=1)
is_msi = is_msi[is_msi["CCLE.MSI.call"].isin(["inferred-MSI","inferred-MSS"])]
is_msi["MSI"] = is_msi["CCLE.MSI.call"] == "inferred-MSI"
```

```python
is_msi.to_hdf("processed/depmap/CCLE_MSI.h5",key="is_msi")
```

# Cell line metadata

```python
# metadata = pd.read_csv("raw/depmap/sample_info.csv",index_col=0)
metadata = pd.read_csv("raw/depmap/sample-info-19q2_v2-achiles-sample-info-full.csv",index_col=0)

metadata["display_disease"] = metadata["disease"].apply(lambda x: x.replace("_"," ").capitalize())
metadata["display_disease"] = metadata["display_disease"].apply(lambda x: "Unknown" if x == " " else x)

metadata.to_csv("processed/depmap/sample_info.csv")
```

# DepMap gene expression

```python
depmap_genex = pd.read_csv("raw/depmap/CCLE_expression.csv",index_col=0)

depmap_genex.to_hdf("processed/depmap/CCLE_expression.h5",key="depmap_genex")
```
