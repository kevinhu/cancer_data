import gzip
import re
import tempfile
from collections import defaultdict

import numpy as np
import pandas as pd

from .. import access

# minimum number of valid samples per
# splicing event for ccle_splicing()
MIN_VALID_COUNT = 100

# minimum standard deviation per splicing
# event for ccle_splicing()
MIN_STDEV = 0.01


class Processors:
    """

    CCLE dataset processors.

    """

    @staticmethod
    def ccle_annotations(raw_path):
        """

        Process CCLE cell line annotations.

        Args:
            raw_path (str): the complete path to the
                            raw downloaded file

        Returns:
            Processed DataFrame

        """

        df = pd.read_csv(raw_path, sep="\t")
        df = df.astype(str)

        return df

    @staticmethod
    def ccle_chromatin(raw_path):
        """

        Process CCLE chromatin profiling estimates.

        Args:
            raw_path (str): the complete path to the
                            raw downloaded file

        Returns:
            Processed DataFrame

        """

        df = pd.read_csv(raw_path, index_col=1)

        df = df.iloc[:, 1:]

        return df

    @staticmethod
    def ccle_translocations_svaba(raw_path):
        """

        Process CCLE SvABA calls.

        Args:
            raw_path (str): the complete path to the
                            raw downloaded file

        Returns:
            Processed DataFrame

        """

        df = pd.read_excel(raw_path)

        ccle_annotations = access.load("ccle_annotations")
        ccle_to_depmap = dict(
            zip(ccle_annotations["CCLE_ID"], ccle_annotations["depMapID"])
        )

        string_cols = [
            "CCLE_name",
            "map_id",
            "bp1",
            "bp2",
            "class",
            "gene1",
            "gene2",
            "site1",
            "site2",
            "fusion",
        ]

        for col in string_cols:
            df[col] = df[col].astype(str)

        df["depmap_id"] = df["CCLE_name"].apply(ccle_to_depmap.get)

        return df

    @staticmethod
    def ccle_rppa_info(raw_path):
        """

        Process CCLE RPPA antibody info.

        Args:
            raw_path (str): the complete path to the
                            raw downloaded file

        Returns:
            Processed DataFrame

        """

        df = pd.read_csv(raw_path)
        df = df.astype(str)

        df["format_id"] = (
            df["Target_Genes"].apply(lambda x: x.replace(" ", "-"))
            + "_"
            + df["Antibody_Name"]
        )

        return df

    @staticmethod
    def ccle_rppa(raw_path):
        """

        Process CCLE RPPA values.

        Args:
            raw_path (str): the complete path to the
                            raw downloaded file

        Returns:
            Processed DataFrame

        """

        df = pd.read_csv(raw_path, index_col=0)

        ccle_rppa_info = access.load("ccle_rppa_info")
        antibody_name_map = dict(
            zip(ccle_rppa_info["Antibody_Name"], ccle_rppa_info["format_id"])
        )

        ccle_annotations = access.load("ccle_annotations")
        ccle_to_depmap = dict(
            zip(ccle_annotations["CCLE_ID"], ccle_annotations["depMapID"])
        )

        df.columns = map(antibody_name_map.get, df.columns)
        df.index = df.index.map(ccle_to_depmap.get)

        return df

    @staticmethod
    def ccle_gene_tpm(raw_path):
        """

        Process CCLE gene TPM measurements.

        Args:
            raw_path (str): the complete path to the
                            raw downloaded file

        Returns:
            Processed DataFrame

        """

        df = pd.read_csv(raw_path, sep="\t", index_col=0)
        df = df.iloc[:, 1:]

        g19_7_definitions = access.load("g19_7_definitions")

        gene_name_map = dict(
            zip(g19_7_definitions["gene_id"], g19_7_definitions["gene_name"])
        )
        gene_name_map = defaultdict(str, gene_name_map)

        gene_names = df.index.map(lambda x: f"{gene_name_map.get(x)}_{x}")

        df = df.set_index(gene_names)
        df = np.log2(df + 1)
        df = df.astype(np.float16)

        df = df.T

        ccle_annotations = access.load("ccle_annotations")
        ccle_to_depmap = dict(
            zip(ccle_annotations["CCLE_ID"], ccle_annotations["depMapID"])
        )

        df.index = df.index.map(ccle_to_depmap.get)

        return df

    @staticmethod
    def ccle_transcript_tpm(raw_path):
        """

        Process CCLE transcript TPM measurements.

        Args:
            raw_path (str): the complete path to the
                            raw downloaded file

        Returns:
            Processed DataFrame

        """

        df = pd.read_csv(raw_path, sep="\t")

        g19_7_definitions = access.load("g19_7_definitions")

        gene_name_map = dict(
            zip(g19_7_definitions["gene_id"], g19_7_definitions["gene_name"])
        )
        gene_name_map = defaultdict(str, gene_name_map)

        df.index = df[["gene_id", "transcript_id"]].apply(
            lambda x: f"{gene_name_map.get(x['gene_id'])}_{x['transcript_id']}", axis=1
        )

        df = df.drop(["gene_id", "transcript_id"], axis=1)

        df = np.log2(df + 1)
        df = df.astype(np.float16)
        df = df.T

        ccle_annotations = access.load("ccle_annotations")
        ccle_to_depmap = dict(
            zip(ccle_annotations["CCLE_ID"], ccle_annotations["depMapID"])
        )

        df.index = df.index.map(ccle_to_depmap.get)

        return df

    @staticmethod
    def ccle_exonusage(raw_path):
        """

        Process CCLE splicing measurements.

        Args:
            raw_path (str): the complete path to the
                            raw downloaded file

        Returns:
            Processed DataFrame

        """

        def reorder_exon(exon):
            """

            Helper function for formatting exon
            splicing IDs.

            Args:
                exon (str): exon identifier

            Returns:
                Processed DataFrame

            """

            exon_split = exon.split("_")
            return "_".join(exon_split[3:]) + "_" + "_".join(exon_split[:3])

        temp = tempfile.NamedTemporaryFile(mode="wb")

        with gzip.open(raw_path, "rb") as f:

            for line in f:

                line = re.sub(b"[^\S\t\n\r]+NA\t", b"nan\t", line)
                line = re.sub(b"[^\S\t\n\r]+NA\n", b"nan\n", line)

                temp.write(line)

        df = pd.read_csv(temp.name, skiprows=2, index_col=0, sep="\t")

        temp.close()

        df.index = df.index.map(reorder_exon)

        df.index = pd.Series(df.index, index=df.index) + "_" + pd.Series(df["gene_id"])
        df = df.iloc[:, 1:]
        df = df.T

        ccle_annotations = access.load("ccle_annotations")
        ccle_to_depmap = dict(
            zip(ccle_annotations["CCLE_ID"], ccle_annotations["depMapID"])
        )

        df.index = df.index.map(ccle_to_depmap.get)

        exonusage_nans = df.isna().sum(axis=0)

        keep_cols = df.columns[exonusage_nans < len(df) - MIN_VALID_COUNT]

        df = df[keep_cols]

        stdevs = df.std(axis=0)
        df = df[df.columns[stdevs >= MIN_STDEV]]

        df = df.astype(np.float16)

        return df

    @staticmethod
    def ccle_mirna(raw_path):
        """

        Process CCLE miRNA measurements.

        Args:
            raw_path (str): the complete path to the
                            raw downloaded file

        Returns:
            Processed DataFrame

        """

        df = pd.read_csv(raw_path, sep="\t", skiprows=2)

        df.index = df["Description"] + "_" + df["Name"].apply(lambda x: x[1:])

        df = df.iloc[:, 2:]
        df = np.log2(df)
        df = df.T
        df = df.astype(np.float16)

        ccle_annotations = access.load("ccle_annotations")
        ccle_to_depmap = dict(
            zip(ccle_annotations["CCLE_ID"], ccle_annotations["depMapID"])
        )

        df.index = df.index.map(ccle_to_depmap.get)

        return df

    @staticmethod
    def ccle_rrbs_tss1kb(raw_path):
        """

        Process CCLE RRBS TSS 1kb measurements.

        Args:
            raw_path (str): the complete path to the
                            raw downloaded file

        Returns:
            Processed DataFrame

        """

        df = pd.read_csv(raw_path, sep="\t", index_col=0)
        df = df.iloc[:-1, 2:]
        df = df.T

        df[df == "\tNA"] = np.nan
        df[df == "    NA"] = np.nan
        df[df == "     NA"] = np.nan
        df = df.astype(np.float16)

        ccle_annotations = access.load("ccle_annotations")
        ccle_to_depmap = dict(
            zip(ccle_annotations["CCLE_ID"], ccle_annotations["depMapID"])
        )

        df.index = df.index.map(ccle_to_depmap.get)

        return df

    @staticmethod
    def ccle_rrbs_tss_clusters(raw_path):
        """

        Process CCLE RRBS TSS cluster measurements.

        Args:
            raw_path (str): the complete path to the
                            raw downloaded file

        Returns:
            Processed DataFrame

        """

        df = pd.read_csv(raw_path, sep="\t", index_col=0)
        df = df.iloc[:-1, 2:]
        df = df.T

        df[df == "\tNA"] = np.nan
        df[df == "    NA"] = np.nan
        df[df == "     NA"] = np.nan
        df = df.astype(np.float16)

        ccle_annotations = access.load("ccle_annotations")
        ccle_to_depmap = dict(
            zip(ccle_annotations["CCLE_ID"], ccle_annotations["depMapID"])
        )

        df.index = df.index.map(ccle_to_depmap.get)

        return df

    @staticmethod
    def ccle_rrbs_cgi_clusters(raw_path):
        """

        Process CCLE RRBS CGI cluster measurements.

        Args:
            raw_path (str): the complete path to the
                            raw downloaded file

        Returns:
            Processed DataFrame

        """

        df = pd.read_csv(raw_path, sep="\t", index_col=0)
        df = df.iloc[:-1]

        df["cluster_pos"] = df.index
        df["cluster_n"] = df.groupby("cluster_pos").cumcount() + 1
        df.index = df["cluster_pos"].astype(str) + "-" + df["cluster_n"].astype(str)

        df = df.iloc[:, 2:-2]
        df = df.T

        df[df == "\tNA"] = np.nan
        df[df == "    NA"] = np.nan
        df[df == "     NA"] = np.nan
        df = df.astype(np.float16)

        ccle_annotations = access.load("ccle_annotations")
        ccle_to_depmap = dict(
            zip(ccle_annotations["CCLE_ID"], ccle_annotations["depMapID"])
        )

        df.index = df.index.map(ccle_to_depmap.get)

        return df

    @staticmethod
    def ccle_rrbs_enhancer_clusters(raw_path):
        """

        Process CCLE RRBS enhancer cluster measurements.

        Args:
            raw_path (str): the complete path to the
                            raw downloaded file

        Returns:
            Processed DataFrame

        """

        df = pd.read_csv(raw_path, sep="\t", index_col=0)

        df.index = df.index + "_" + df.groupby(level=0).cumcount().astype(str)

        df = df.iloc[:, 2:]
        df.index = df.index.map(lambda x: x.replace("_", "-")) + "_enh"
        df = df.T

        df[df == "\tNA"] = np.nan
        df[df == "    NA"] = np.nan
        df[df == "     NA"] = np.nan
        df = df.astype(np.float16)

        ccle_annotations = access.load("ccle_annotations")
        ccle_to_depmap = dict(
            zip(ccle_annotations["CCLE_ID"], ccle_annotations["depMapID"])
        )

        df.index = df.index.map(ccle_to_depmap.get)

        return df

    @staticmethod
    def ccle_tertp(raw_path):
        """

        Process CCLE TERT promoter calls.

        Args:
            raw_path (str): the complete path to the
                            raw downloaded file

        Returns:
            Processed DataFrame

        """

        df = pd.read_excel(raw_path, skiprows=4)

        df = df.set_index("depMapID")
        df["TERTp_mut"] = df["TERT_promoter_mutation"] != "wildtype"

        return df

    @staticmethod
    def ccle_msi(raw_path):
        """

        Process CCLE MSI calls.

        Args:
            raw_path (str): the complete path to the
                            raw downloaded file

        Returns:
            Processed DataFrame

        """

        df = pd.read_excel(raw_path, sheet_name="MSI calls")

        df = df[df["CCLE.MSI.call"].isin(["inferred-MSI", "inferred-MSS"])]

        df = df.astype(str)

        df["MSI"] = df["CCLE.MSI.call"] == "inferred-MSI"

        df = df.set_index("depMapID")

        return df

    @staticmethod
    def ccle_metabolomics(raw_path):
        """

        Process CCLE metabolomics measurements.

        Args:
            raw_path (str): the complete path to the
                            raw downloaded file

        Returns:
            Processed DataFrame

        """

        df = pd.read_csv(raw_path)

        df["DepMap_ID"] = df["DepMap_ID"].astype(str)

        df = df.set_index("DepMap_ID")
        df = df.drop(["CCLE_ID"], axis=1)

        df = df.astype(np.float16)

        return df

    @staticmethod
    def ccle_proteomics(raw_path):
        """

        Process CCLE proteomics measurements.

        Args:
            raw_path (str): the complete path to the
                            raw downloaded file

        Returns:
            Processed DataFrame

        """

        df = pd.read_csv(raw_path)
        df.index = df["Gene_Symbol"].fillna("UNNAMED") + "_" + df["Uniprot_Acc"]

        df = df.iloc[:, 48:]

        df = df.T

        ccle_annotations = access.load("ccle_annotations")
        ccle_to_depmap = dict(
            zip(ccle_annotations["CCLE_ID"], ccle_annotations["depMapID"])
        )

        df.index = df.index.map(lambda x: "_".join(x.split("_")[:-1]))
        df.index = df.index.map(ccle_to_depmap.get)

        df = df[~df.index.duplicated(keep="first")]

        df = df.astype(np.float16)

        return df
