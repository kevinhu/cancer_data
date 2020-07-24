# cancer-data

This package provides that handle the downloading and processing of several public genomics datasets useful for cancer research.

## Datasets

A complete description of the datasets may be found in the [schema](https://github.com/kevinhu/cancer-data/blob/master/cancer_data/schema.txt).

| Collection                                    | Datasets                                                     | Portal                                                       |
| --------------------------------------------- | ------------------------------------------------------------ | ------------------------------------------------------------ |
| Cancer Cell Line Encyclopedia (CCLE)          | Many (see portal)                                            | https://portals.broadinstitute.org/ccle/data (registration required) |
| Cancer Dependency Map (DepMap)                | Genome-wide CRISPR-cas9 and RNAi screens, gene expression, mutations, and copy number | https://depmap.org/portal/download/                          |
| The Cancer Genome Atlas (TCGA)                | Mutations, RNAseq expression and splicing, and copy number   | https://xenabrowser.net/datapages/?cohort=TCGA%20Pan-Cancer%20(PANCAN)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443 |
| The Genotype-Tissue Expression (GTEx) Project | RNAseq expression and splicing                               | https://gtexportal.org/home/datasets                         |

## Features

The goal of this package is to make statistical analysis and coordination of these datasets easier. To that end, it provides the following features:

1. Harmonization: datasets within a collection have sample IDs reduced to the same format. For instance, all CCLE+DepMap datasets have been modified to use Achilles/Arxspan IDs, rather than cell line names.
2. Speed: processed datasets are all stored in high-performance [HDF5 format](https://en.wikipedia.org/wiki/Hierarchical_Data_Format), allowing large tables to be loaded magnitudes faster than in raw CSV or TSV formats.
3. Space: when possible, tables of purely numerical values (e.g. gene expression, methylation, drug sensitivities) are stored in half-precision format.