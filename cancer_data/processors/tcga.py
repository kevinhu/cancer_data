import numpy as np
import pandas as pd

from ..utils import concat_cols

# minimum number of valid samples per
# splicing event for tcga_splicing()
MIN_VALID_COUNT = 100

# minimum standard deviation per splicing
# event for tcga_splicing()
MIN_STDEV = 0.01

# tumor abbreviations map
TCGA_MAP = {
    "acute myeloid leukemia": "LAML",
    "adrenocortical cancer": "ACC",
    "bladder urothelial carcinoma": "BLCA",
    "brain lower grade glioma": "LGG",
    "breast invasive carcinoma": "BRCA",
    "cervical & endocervical cancer": "CESC",
    "cholangiocarcinoma": "CHOL",
    "colon adenocarcinoma": "COAD",
    "diffuse large B-cell lymphoma": "DLBC",
    "esophageal carcinoma": "ESCA",
    "glioblastoma multiforme": "GBM",
    "head & neck squamous cell carcinoma": "HNSC",
    "kidney chromophobe": "KICH",
    "kidney clear cell carcinoma": "KIRC",
    "kidney papillary cell carcinoma": "KIRP",
    "liver hepatocellular carcinoma": "LIHC",
    "lung adenocarcinoma": "LUAD",
    "lung squamous cell carcinoma": "LUSC",
    "mesothelioma": "MESO",
    "ovarian serous cystadenocarcinoma": "OV",
    "pancreatic adenocarcinoma": "PAAD",
    "pheochromocytoma & paraganglioma": "PCPG",
    "prostate adenocarcinoma": "PRAD",
    "rectum adenocarcinoma": "READ",
    "sarcoma": "SARC",
    "skin cutaneous melanoma": "SKCM",
    "stomach adenocarcinoma": "STAD",
    "testicular germ cell tumor": "TGCT",
    "thymoma": "THYM",
    "thyroid carcinoma": "THCA",
    "uterine carcinosarcoma": "UCS",
    "uterine corpus endometrioid carcinoma": "UCEC",
    "uveal melanoma": "UVM",
}


def tcga_splicing(raw_path):
    """

    General handler for all TCGA splicing files.

    Args:
        raw_path (str): the complete path to the
                        raw downloaded file

    Returns:
        Processed DataFrame

    """

    chunk_iterator = pd.read_csv(raw_path, sep="\t", chunksize=1000)

    merged = []

    for chunk in chunk_iterator:

        chunk["exon_id"] = concat_cols(
            chunk,
            [
                "gene_name",
                "event_type",
                "event_chr",
                "event_coordinates",
                "alt_region_coordinates",
            ],
            "_",
        )

        chunk = chunk.drop(
            [
                "event_id",
                "event_type",
                "event_chr",
                "event_coordinates",
                "alt_region_coordinates",
                "gene_name",
            ],
            axis=1,
        )

        chunk = chunk.set_index("exon_id")
        chunk = chunk.astype(np.float16)
        chunk = chunk.T

        nan_counts = chunk.isna().sum(axis=0)

        keep_cols = chunk.columns[nan_counts < len(chunk) - MIN_VALID_COUNT]

        chunk = chunk.filter(keep_cols, axis=1)

        stdevs = chunk.std(axis=0)

        keep_cols = chunk.columns[stdevs >= MIN_STDEV]

        chunk = chunk.filter(keep_cols, axis=1)

        merged.append(chunk)

    merged = pd.concat(merged, axis=1)

    # remove prefix identifiers from names
    merged.index = merged.index.map(lambda x: x.split(".")[0])

    return merged


class Processors:
    """

    TCGA dataset processors.

    """

    @staticmethod
    def tcga_annotations(raw_path):
        """

        Process TCGA sample annotations.

        Args:
            raw_path (str): the complete path to the
                            raw downloaded file

        Returns:
            Processed DataFrame

        """

        df = pd.read_csv(raw_path, sep="\t", index_col=0)

        df["abbreviated_disease"] = df["_primary_disease"].apply(lambda x: TCGA_MAP[x])

        str_cols = ["sample_type", "_primary_disease", "abbreviated_disease"]

        for col in str_cols:

            df[col] = df[col].astype(str)

        df.index = df.index.astype(str)

        return df

    @staticmethod
    def tcga_mutations(raw_path):
        """

        Process TCGA MC3 mutation calls.

        Args:
            raw_path (str): the complete path to the
                            raw downloaded file

        Returns:
            Processed DataFrame

        """

        df = pd.read_csv(raw_path, sep="\t")

        df = df.astype(str)
        df["start"] = df["start"].astype(int)
        df["end"] = df["end"].astype(int)

        df["DNA_VAF"] = df["DNA_VAF"].astype(np.float16)

        return df

    @staticmethod
    def tcga_msi(raw_path):
        """

        Process TCGA MSI calls.

        Args:
            raw_path (str): the complete path to the
                            raw downloaded file

        Returns:
            Processed DataFrame

        """

        df = pd.read_excel(raw_path)
        df = df[df["Case ID"].apply(lambda x: x[:4] == "TCGA")]
        df["sample_type"] = df["Tumor Filename"].apply(lambda x: x.split("-")[3][:-1])
        df["Case ID"] = df["Case ID"] + "-" + df["sample_type"]
        df = df.set_index("Case ID")

        return df

    @staticmethod
    def tcga_cn_continuous(raw_path):
        """

        Process TCGA continuous copy number calls.

        Args:
            raw_path (str): the complete path to the
                            raw downloaded file

        Returns:
            Processed DataFrame

        """

        df = pd.read_csv(raw_path, sep="\t", index_col=0)
        df = df.T
        df = df.astype(np.float16)

        return df

    @staticmethod
    def tcga_cn_thresholded(raw_path):
        """

        Process TCGA discrete (thresholded) copy number calls.

        Args:
            raw_path (str): the complete path to the
                            raw downloaded file

        Returns:
            Processed DataFrame

        """

        df = pd.read_csv(raw_path, sep="\t", index_col=0)
        df = df.T
        df = df.astype(np.float16)

        return df

    @staticmethod
    def tcga_cn_whitelisted(raw_path):
        """

        Process TCGA whitelisted copy number calls.

        Args:
            raw_path (str): the complete path to the
                            raw downloaded file

        Returns:
            Processed DataFrame

        """

        df = pd.read_csv(raw_path, sep="\t", index_col=0)
        df = df.T
        df = df.astype(np.float16)

        return df

    @staticmethod
    def tcga_normalized_gene_expression(raw_path):
        """

        Process TCGA batch effects-normalized gene expression
        estimates.

        Args:
            raw_path (str): the complete path to the
                            raw downloaded file

        Returns:
            Processed DataFrame

        """

        df = pd.read_csv(raw_path, sep="\t", index_col=0)

        df = df.T
        df.columns = [df.columns[i] + "_" + str(i) for i in range(len(df.columns))]

        df = df.astype(np.float16)

        return df

    @staticmethod
    def tcga_a3ss(raw_path):
        """

        Process TCGA A3SS splicing events.

        Args:
            raw_path (str): the complete path to the
                            raw downloaded file

        Returns:
            Processed DataFrame

        """

        return tcga_splicing(raw_path)

    @staticmethod
    def tcga_a5ss(raw_path):
        """

        Process TCGA A5SS splicing events.

        Args:
            raw_path (str): the complete path to the
                            raw downloaded file

        Returns:
            Processed DataFrame

        """

        return tcga_splicing(raw_path)

    @staticmethod
    def tcga_se(raw_path):
        """

        Process TCGA SE splicing events.

        Args:
            raw_path (str): the complete path to the
                            raw downloaded file

        Returns:
            Processed DataFrame

        """

        return tcga_splicing(raw_path)

    @staticmethod
    def tcga_ir(raw_path):
        """

        Process TCGA IR splicing events.

        Args:
            raw_path (str): the complete path to the
                            raw downloaded file

        Returns:
            Processed DataFrame

        """

        return tcga_splicing(raw_path)

    @staticmethod
    def tcga_mx(raw_path):
        """

        Process TCGA MX splicing events.

        Args:
            raw_path (str): the complete path to the
                            raw downloaded file

        Returns:
            Processed DataFrame

        """

        return tcga_splicing(raw_path)
