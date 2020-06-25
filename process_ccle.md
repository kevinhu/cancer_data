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

from scipy.stats import zscore

from collections import Counter, defaultdict

import sys
import os
import re
sys.path.append(os.path.relpath("../helper"))

import name_mappings

from gtfparse import read_gtf
```

```python
cell_line_info = pd.read_csv("raw/depmap/DepMap-2018q4-celllines.csv",index_col=0)
```

# Gencode v19 definitions

```python
g19_definitions = read_gtf("raw/ccle/gencode.v19.genes.v7_model.patched_contigs.gtf")

g19_definitions.to_hdf("processed/ccle/gencode.v19.genes.v7_model.patched_contigs.h5",
                       key="g19_definitions",
                       mode="w")
```

```python
g19_gene_definitions = g19_definitions[g19_definitions["feature"]=="gene"]

ensembl_id_map = dict(zip(g19_gene_definitions["gene_id"],g19_gene_definitions["gene_name"]))
ensembl_id_map = defaultdict(str, ensembl_id_map)
```

# ENSEMBL 75 definitions

```python
ensembl_definitions = read_gtf("raw/ccle/Homo_sapiens.GRCh37.75.gtf")

ensembl_definitions.to_hdf("processed/ccle/Homo_sapiens.GRCh37.75.h5",
                       key="ensembl_definitions",
                       mode="w")
```

# Chromatin profiling

```python
chromatin_profiling = pd.read_csv("raw/ccle/CCLE_GlobalChromatinProfiling_20181130.csv",index_col=1).iloc[:,1:]
chromatin_profiling.to_hdf("processed/ccle/CCLE_GlobalChromatinProfiling_20181130.hdf",
                           key="chromatin_profiling",
                           mode="w")

```

# RPPA

```python
rppa = pd.read_csv("raw/ccle/CCLE_RPPA_20181003.csv",index_col=0)
rppa.index = rppa.index.map(lambda x: name_mappings.name_map[x])
rppa.to_hdf("processed/ccle/CCLE_RPPA_20181003.hdf",
            key="rppa",
            mode="w")
```

# RNAseq


## Exonusage

```python
def reorder_exon(exon):
    exon_split = exon.split("_")
    return "_".join(exon_split[3:]) + "_" + "_".join(exon_split[:3])

exonusage = pd.read_csv("raw/ccle/CCLE_RNAseq_ExonUsageRatio_20180929.gct",skiprows=2,index_col=0,sep="\t")
exonusage.index = exonusage.index.map(lambda x:reorder_exon(x))
```

```python
exonusage[exonusage=="\tNA"] = np.nan
exonusage[exonusage=="    NA"] = np.nan
exonusage[exonusage=="     NA"] = np.nan
```

```python
exon_ids = pd.Series(exonusage.index,index=exonusage.index) + "_" + pd.Series(exonusage["gene_id"])
exonusage = exonusage.set_index(exon_ids).iloc[:,1:]
exonusage = exonusage.T
exonusage.index = exonusage.index.map(lambda x: name_mappings.name_map[x])
```

```python
exonusage = exonusage.astype(np.float32)
exonusage_nans = exonusage.isna().sum(axis=0)
sns.distplot(exonusage_nans)
exonusage = exonusage[exonusage.columns[exonusage_nans<len(exonusage)-100]]
```

```python
exonusage_stdevs = exonusage.std(axis=0)
sns.distplot(exonusage_stdevs)
exonusage = exonusage[exonusage.columns[exonusage_stdevs>0.01]]
```

```python
exonusage.to_hdf("processed/ccle/CCLE_RNAseq_ExonUsageRatio_20180929.hdf",
                 key="exonusage",
                 mode="w")
```

## Intron retention

```python
intronret = pd.read_csv("raw/ccle/CCLE_intron_retention.csv",index_col=0)
intronret.columns = ["_".join(x.split(";")[::-1]) for x in intronret.columns]
intronret.index = intronret.index.map(lambda x: name_mappings.name_map[x])
```

```python
intronret_nans = intronret.isna().sum(axis=0)
sns.distplot(intronret_nans)
intronret = intronret[intronret.columns[intronret_nans<800]]
```

```python
intronret_stdevs = intronret.std(axis=0)
sns.distplot(intronret_stdevs)
intronret = intronret[intronret.columns[intronret_stdevs>0.1]]
```

```python
intronret.to_hdf("processed/ccle/CCLE_intron_retention.hdf",
                 key="intronret",
                 mode="w")
```

## Transcript expression

```python
ccle_transcripts = pd.read_csv("raw/ccle/CCLE_RNAseq_rsem_transcripts_tpm_20180929.txt",sep="\t",index_col=1)
ccle_transcripts["gene_id"] = ccle_transcripts["gene_id"].apply(lambda x: ensembl_id_map[x])

```

```python
gene_transcript_ids = ccle_transcripts["gene_id"] + "_" + pd.Series(ccle_transcripts.index,index=ccle_transcripts.index)
ccle_transcripts = ccle_transcripts.set_index(gene_transcript_ids)
```

```python
ccle_transcripts = ccle_transcripts.iloc[:,1:]
ccle_transcripts = np.log2(ccle_transcripts+1)
ccle_transcripts = ccle_transcripts.T
ccle_transcripts.index = ccle_transcripts.index.map(lambda x: name_mappings.name_map[x])
```

```python
ccle_transcript_stdevs = ccle_transcripts.std(axis=0)
sns.distplot(ccle_transcript_stdevs)
ccle_transcripts = ccle_transcripts[ccle_transcripts.columns[ccle_transcript_stdevs>0.25]]
```

```python
ccle_transcripts.to_hdf("processed/ccle/CCLE_RNAseq_rsem_transcripts_tpm_20180929.hdf",
                        key="ccle_transcripts",
                        mode="w")
```

## Gene expression

```python
ccle_genex = pd.read_csv("raw/ccle/CCLE_RNAseq_rsem_genes_tpm_20180929.txt",sep="\t",index_col=0)

ccle_genex = ccle_genex.iloc[:,1:]
ccle_gene_names = ccle_genex.index.map(lambda x: ensembl_id_map[x])
gene_names_ids = ccle_gene_names + "_" + pd.Series(ccle_genex.index,index=ccle_genex.index)
ccle_genex = ccle_genex.set_index(gene_names_ids)
ccle_genex = np.log2(ccle_genex+1)
ccle_genex = ccle_genex.T
ccle_genex.index = ccle_genex.index.map(lambda x: name_mappings.name_map[x])
```

```python
ccle_genex.to_hdf("processed/ccle/CCLE_RNAseq_rsem_genes_tpm_20180929.hdf",
                  key="ccle_genex",
                  mode="w")
```

## RPKM expression

```python
ccle_genex = pd.read_csv("raw/ccle/CCLE_RNAseq_genes_rpkm_20180929.gct",sep="\t",index_col=0,skiprows=2)

gene_names_ids = ccle_genex["Description"] + "_" + ccle_genex.index
ccle_genex = ccle_genex.iloc[:,1:]
ccle_genex = ccle_genex.set_index(gene_names_ids)
ccle_genex = np.log2(ccle_genex+1)
ccle_genex = ccle_genex.T
ccle_genex.index = ccle_genex.index.map(lambda x: name_mappings.name_map[x])
```

```python
ccle_genex.to_hdf("processed/ccle/CCLE_RNAseq_genes_rpkm_20180929.hdf",
                  key="ccle_genex",
                  mode="w")
```

## ssGSEA

```python
name_map = dict(zip(cell_line_info["CCLE_Name"],cell_line_info.index))
```

```python
ssgsea = pd.read_csv("raw/ccle/GSEA_RNAseq.csv",index_col=0)
ssgsea = ssgsea[ssgsea.index.isin(name_map.keys())]
ssgsea.index = ssgsea.index.map(lambda x: name_map[x])
ssgsea = ssgsea.dropna(thresh=50,axis=0)
```

```python
ssgsea.to_hdf("processed/ccle/GSEA_RNAseq.hdf",
              key="ssgsea",
              mode="w")
```

# RRBS


## CpG islands

```python
cgi_meth = pd.read_csv("raw/ccle/CCLE_RRBS_cgi_CpG_clusters_20181119.txt",sep="\t",index_col=0)
cgi_meth = cgi_meth.iloc[:-1]
cgi_meth["cluster_pos"] = cgi_meth.index
cgi_meth['cluster_n'] = cgi_meth.groupby('cluster_pos').cumcount()+1
cgi_meth.index = cgi_meth["cluster_pos"].astype(str) + "-" + cgi_meth["cluster_n"].astype(str)

cgi_meth = cgi_meth.iloc[:,2:-2]
cgi_meth = cgi_meth.T
cgi_meth = cgi_meth.astype(np.float32)
cgi_meth.index = cgi_meth.index.map(lambda x: name_mappings.name_map[x])
```

```python
cgi_meth_stds = cgi_meth.std(axis=0)
sns.distplot(cgi_meth_stds)
cgi_meth = cgi_meth[cgi_meth.columns[cgi_meth_stds>0.05]]
```

```python
cgi_meth.to_hdf("processed/ccle/CCLE_RRBS_cgi_CpG_clusters_20181119.hdf",
                key="cgi_meth",
                mode="w")
```

## Enhancers

```python
enh_meth = pd.read_csv("raw/ccle/CCLE_RRBS_enh_CpG_clusters_20181119.txt",sep="\t",index_col=0)
enh_meth.index = enh_meth.index + "_" + enh_meth.groupby(level=0).cumcount().astype(str)

enh_meth = enh_meth.iloc[:,2:]
enh_meth.index = enh_meth.index.map(lambda x: x.replace("_","-"))+"_enh"
enh_meth = enh_meth.T
enh_meth.index = enh_meth.index.map(lambda x: name_mappings.name_map[x])
```

```python
enh_meth[enh_meth=="\tNA"] = np.nan
enh_meth[enh_meth=="    NA"] = np.nan
enh_meth[enh_meth=="     NA"] = np.nan
enh_meth = enh_meth.astype(float)
```

```python
enh_meth_stds = enh_meth.std(axis=0)
sns.distplot(enh_meth_stds)
enh_meth = enh_meth[enh_meth.columns[enh_meth_stds>0.05]]
```

```python
enh_meth.to_hdf("processed/ccle/CCLE_RRBS_enh_CpG_clusters_20181119.hdf",
                key="enh_meth",
                mode="w")
```

## TSS 1kb

```python
tss1kb_meth = pd.read_csv("raw/ccle/CCLE_RRBS_TSS1kb_20181022.txt",sep="\t",index_col=0)
tss1kb_meth = tss1kb_meth.iloc[:-1,2:]
tss1kb_meth = tss1kb_meth.T
tss1kb_meth.index = tss1kb_meth.index.map(lambda x: name_mappings.name_map[x])
```

```python
tss1kb_meth[tss1kb_meth=="\tNA"] = np.nan
tss1kb_meth[tss1kb_meth=="    NA"] = np.nan
tss1kb_meth[tss1kb_meth=="     NA"] = np.nan
tss1kb_meth = tss1kb_meth.astype(float)
```

```python
tss1kb_meth_stds = tss1kb_meth.std(axis=0)
sns.distplot(tss1kb_meth_stds)
tss1kb_meth = tss1kb_meth[tss1kb_meth.columns[tss1kb_meth_stds>0.05]]
```

```python
tss1kb_meth.to_hdf("processed/ccle/CCLE_RRBS_TSS1kb_20181022.hdf",
                   key="tss1kb_meth",
                   mode="w"
                  )
```

## TSS clusters

```python
tssclust_meth = pd.read_csv("raw/ccle/CCLE_RRBS_tss_CpG_clusters_20181022.txt",sep="\t",index_col=0)
tssclust_meth = tssclust_meth.iloc[:-1,2:]
tssclust_meth = tssclust_meth.T
tssclust_meth.index = tssclust_meth.index.map(lambda x: name_mappings.name_map[x])
```

```python
tssclust_meth[tssclust_meth=="\tNA"] = np.nan
tssclust_meth[tssclust_meth=="    NA"] = np.nan
tssclust_meth[tssclust_meth=="     NA"] = np.nan
tssclust_meth = tssclust_meth.astype(float)
```

```python
tssclust_meth_stds = tssclust_meth.std(axis=0)
sns.distplot(tssclust_meth_stds)
tssclust_meth = tssclust_meth[tssclust_meth.columns[tssclust_meth_stds>0.05]]
```

```python
tssclust_meth.to_hdf("processed/ccle/CCLE_RRBS_tss_CpG_clusters_20181022.hdf",key="tssclust_meth")
```

# miRNA

```python
mirna = pd.read_csv("raw/ccle/CCLE_miRNA_20181103.gct.txt",sep="\t",skiprows=2)
mirna.index = mirna["Description"] + "_" + mirna["Name"].apply(lambda x: x[1:])

mirna = mirna.iloc[:,2:]
mirna = np.log2(mirna.T)

mirna.index = mirna.index.map(lambda x: name_mappings.name_map[x])
```

```python
mirna.to_hdf("processed/ccle/CCLE_miRNA_20181103.hdf", key="mirna", mode="w")
```

# TERT promoter

```python
tertp_info = pd.read_csv("raw/ccle/Supplementary_Table_2_TERT_promoter_mutation.csv",index_col=0)

tertp_info.index = tertp_info.index.map(lambda x: name_mappings.name_map[x])

tertp_info["tertp_mut"] = tertp_info["Classical_Promoter_Mutations"] != "wildtype"

tertp_info.to_hdf("processed/ccle/Supplementary_Table_2_TERT_promoter_mutation.h5",
                  key="tertp_info",
                  mode="w")
```

# WES/WGS mutation calls

```python
genome_muts = pd.read_csv("raw/ccle/ccle2maf_ExcludSangerDriftedSubset_20180820.txt",
                          sep="\t")
```

# Drug sensitivities

```python
primary_mad = pd.read_csv("raw/ccle/primary-screen-public-tentative_v4-primary-merged-mad-lfcvc-cb.csv",
                         index_col=0)
primary_median = pd.read_csv("raw/ccle/primary-screen-public-tentative_v4-primary-merged-median-lfcvc-cb.csv",
                            index_col=0)
primary_meta = pd.read_csv("raw/ccle/primary-screen-public-tentative_v4-primary-merged-row-meta.csv")

secondary_mad = pd.read_csv("raw/ccle/secondary-screen-public-tentative_v3-secondary-merged-mad-lfcvc-cb.csv",
                           index_col=0)
secondary_median = pd.read_csv("raw/ccle/secondary-screen-public-tentative_v3-secondary-merged-median-lfcvc-cb.csv",
                              index_col=0)
secondary_meta = pd.read_csv("raw/ccle/secondary-screen-public-tentative_v3-secondary-merged-row-meta.csv")

primary_name_map = dict(zip(primary_meta["feature_id"],primary_meta["depmap_id"]))
secondary_name_map = dict(zip(secondary_meta["feature_id"],secondary_meta["depmap_id"]))


```

```python
primary_mad.index = primary_mad.index.map(lambda x: primary_name_map[x])
primary_median.index = primary_median.index.map(lambda x: primary_name_map[x])

secondary_mad.index = secondary_mad.index.map(lambda x: secondary_name_map[x])
secondary_median.index = secondary_median.index.map(lambda x: secondary_name_map[x])
```

```python
primary_mad.to_hdf(
    "processed/ccle/primary-screen-public-tentative_v4-primary-merged-mad-lfcvc-cb.h5", 
    key="primary_mad",
    mode="w")
primary_median.to_hdf(
    "processed/ccle/primary-screen-public-tentative_v4-primary-merged-median-lfcvc-cb.h5", 
    key="primary_median",
    mode="w")

secondary_mad.to_hdf(
    "processed/ccle/secondary-screen-public-tentative_v3-secondary-merged-mad-lfcvc-cb.h5", 
    key="secondary_mad",
    mode="w")
secondary_median.to_hdf(
    "processed/ccle/secondary-screen-public-tentative_v3-secondary-merged-median-lfcvc-cb.h5", 
    key="secondary_median",
    mode="w")
```

# Raw mutation data

```python
mutation_calls = pd.read_csv("../data/raw/depmap/CCLE_mutations_19q4.csv",index_col=0)
```

```python
damaging_muts = mutation_calls[mutation_calls["Variant_annotation"]=="damaging"]
hs_muts = mutation_calls[(mutation_calls["isCOSMIChotspot"]==True)|(mutation_calls["isTCGAhotspot"]==True)]

damaging_counts = Counter(damaging_muts["Hugo_Symbol"])
hs_counts = Counter(hs_muts["Hugo_Symbol"])

damaging_muts["count"] = damaging_muts["Hugo_Symbol"].apply(lambda x: damaging_counts[x])
hs_muts["count"] = hs_muts["Hugo_Symbol"].apply(lambda x: hs_counts[x])

damaging_muts = damaging_muts[damaging_muts["count"]>=4]
hs_muts = hs_muts[hs_muts["count"]>=4]

damaging_muts["id"] = damaging_muts["Hugo_Symbol"] + "_" + damaging_muts["DepMap_ID"]
hs_muts["id"] = hs_muts["Hugo_Symbol"] + "_" + hs_muts["DepMap_ID"]

damaging_muts = damaging_muts.drop_duplicates(subset=["id"],keep="first")
hs_muts = hs_muts.drop_duplicates(subset=["id"],keep="first")
```

```python
hs_muts["value"] = 1
damaging_muts["value"] = 1

hs_mut_mat = pd.pivot_table(hs_muts, values="value", index=[
                            "DepMap_ID"], columns="Hugo_Symbol", fill_value=0)
damaging_mut_mat = pd.pivot_table(damaging_muts, values="value", index=[
                                  "DepMap_ID"], columns="Hugo_Symbol", fill_value=0)
```

```python
hs_mut_mat.to_hdf("processed/ccle/hs_mut_mat.hdf",key="hs_mut_mat",mode="w")
damaging_mut_mat.to_hdf("processed/ccle/damaging_mut_mat.hdf",key="damaging_mut_mat",mode="w")
```

# Binary calls for CNA/mutations

```python
revealer = pd.read_csv("raw/ccle/CCLE_MUT_CNA_AMP_DEL_binary_Revealer.gct",
                       skiprows=2, index_col=0, sep="\t").iloc[:, 1:].T
revealer = revealer.drop("UMRC6_KIDNEY",axis=0)

revealer.index = revealer.index.map(lambda x: name_mappings.name_map[x])
revealer.to_hdf("processed/ccle/CCLE_MUT_CNA_AMP_DEL_binary_Revealer.hdf",key="revealer",mode="w")
```

# ssGSEA (new)

```python
oncogenic = pd.read_csv("raw/ccle/CCLE_RNAseq_genes_rpkm_20180929.ssgsea_ONCOGENIC.gct",
                       sep="\t",
                       index_col=0,
                       skiprows=2
                      )
hallmark = pd.read_csv("raw/ccle/CCLE_RNAseq_genes_rpkm_20180929.ssgsea_HALLMARK.gct",
                       sep="\t",
                       index_col=0,
                       skiprows=2
                      )
senescence = pd.read_csv("raw/ccle/CCLE_RNAseq_genes_rpkm_20180929.ssgsea_SENESCENCE.gct",
                       sep="\t",
                       index_col=0,
                       skiprows=2
                      )
oncogenic = oncogenic.iloc[:,1:].T
hallmark = hallmark.iloc[:,1:].T
senescence = senescence.iloc[:,1:].T

oncogenic = oncogenic.apply(zscore)
hallmark = hallmark.apply(zscore)
senescence = senescence.apply(zscore)

ssgsea = pd.concat([oncogenic,hallmark,senescence],axis=1,sort=True)
ssgsea.index = ssgsea.index.map(lambda x: name_mappings.name_map[x])
ssgsea.to_hdf("processed/ccle/ssgsea.hdf",key="ssgsea",mode="w")
```

# Protein levels

```python
ms_prot = pd.read_excel("raw/ccle/Table_S2_Protein_Quant_Normalized.xlsx",
                        sheet_name="Normalized Protein Expression")
ms_prot = ms_prot[[x for x in ms_prot if x[:6] != "Column"]]
ms_prot["Protein_Id"] = ms_prot["Uniprot"] + "_" + ms_prot["Uniprot_Acc"]
ms_prot = ms_prot.set_index("Protein_Id")
ms_prot = ms_prot.iloc[:, 47:]
ms_prot = ms_prot.T

ms_prot.index = ms_prot.index.map(lambda x: "_".join(x.split("_")[:-1]))
ms_prot.index = ms_prot.index.map(lambda x: name_mappings.name_map[x])

ms_prot = ms_prot[~ms_prot.index.duplicated(keep="first")]
```

```python
protein_mapping = pd.read_csv("raw/ccle/HUMAN_9606_idmapping.dat",sep="\t",names=["uniprot_id","feature","converted_id"])


```

```python
ensembl_protein_gene_mapping = protein_mapping[protein_mapping["feature"]=="Ensembl"]
```

```python
ensembl_protein_gene_mapping.shape
```

```python
len(set(ensembl_protein_gene_mapping["uniprot_id"]))
```

```python
ms_prot.to_hdf("../data/processed/ccle/ms_prot.h5", key="ms_prot", mode="w")
```
