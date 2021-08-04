# cancer_data

This package provides unified methods for accessing popular datasets used in cancer research.

**[Full documentation](https://cancer_data.kevinhu.io)**

## Installation

```bash
pip install cancer_data
```

## System requirements

The raw downloaded files occupy approximately 15 GB, and the processed HDFs take up about 10 GB. On a relatively recent machine with a fast SSD, processing all of the files after download takes about 3-4 hours. At least 16 GB of RAM is recommended for handling the large splicing tables.

## Datasets

A complete description of the datasets may be found in [schema.csv](https://github.com/kevinhu/cancer-data/blob/master/cancer_data/schema.csv).

| Collection                                    | Datasets                                                                              | Portal                                                                                                                          |
| --------------------------------------------- | ------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------- |
| Cancer Cell Line Encyclopedia (CCLE)          | Many (see portal)                                                                     | https://portals.broadinstitute.org/ccle/data (registration required)                                                            |
| Cancer Dependency Map (DepMap)                | Genome-wide CRISPR-cas9 and RNAi screens, gene expression, mutations, and copy number | https://depmap.org/portal/download/                                                                                             |
| The Cancer Genome Atlas (TCGA)                | Mutations, RNAseq expression and splicing, and copy number                            | https://xenabrowser.net/datapages/?cohort=TCGA%20Pan-Cancer%20(PANCAN)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443 |
| The Genotype-Tissue Expression (GTEx) Project | RNAseq expression and splicing                                                        | https://gtexportal.org/home/datasets                                                                                            |

## Features

The goal of this package is to make statistical analysis and coordination of these datasets easier. To that end, it provides the following features:

1. Harmonization: datasets within a collection have sample IDs reduced to the same format. For instance, all CCLE+DepMap datasets have been modified to use Achilles/Arxspan IDs, rather than cell line names.
2. Speed: processed datasets are all stored in high-performance [HDF5 format](https://en.wikipedia.org/wiki/Hierarchical_Data_Format), allowing large tables to be loaded orders of magnitude faster than with CSV or TSV formats.
3. Space: tables of purely numerical values (e.g. gene expression, methylation, drug sensitivities) are stored in half-precision format. Compression is used for all tables, resulting in size reductions by factors of over 10 for sparse matrices such as mutation tables, and over 50 for highly-redundant tables such as gene-level copy number estimates.

## How it works

The [schema](https://github.com/kevinhu/cancer-data/blob/master/cancer_data/schema.csv) serves as the reference point for all datasets used. Each dataset is identified by a unique `id` column, which also serves as its access identifier.

Datasets are downloaded from the location specified in `download_url`, after which they are checked against the provided `downloaded_md5` hash.

The next steps depend on the `type` of the dataset:

-   `reference` datasets, such as the hg19 FASTA files, are left as-is.
-   `primary_dataset` objects are preprocessed and converted into HDF5 format.
-   `secondary_dataset` objects are defined as being made from `primary_dataset` objects. These are also processed and converted into HDF5 format.

To keep track of which datasets are necessary for producing another, the `dependencies` column specifies the dataset `id`s that are required for making another. For instance, the `ccle_proteomics` dataset is dependent on the `ccle_annotations` dataset for converting cell line names to Achilles IDs. When running the processing pipeline, the package will automatically check that dependencies are met, and raise an error if they are not found.

## Notes

Some datasets have filtering applied to reduce their size. These are listed below:

-   CCLE, GTEx, and TCGA splicing datasets have been filtered to remove splicing events with many missing values as well as those with low standard deviations.
-   When constructing binary mutation matrices (`depmap_damaging` and `depmap_hotspot`), a minimum mutation frequency is used to remove especially rare (present in less than four samples) mutations.
-   The TCGA MX splicing dataset is extremely large (approximately 10,000 rows by 900,000 columns), so it has been split column-wise into 8 chunks.
