import numpy as np
import pandas as pd

from .. import access
from ..utils import concat_cols

# minimum number of valid samples per
# splicing event for gtex_splicing()
MIN_VALID_COUNT = 100

# minimum standard deviation per splicing
# event for gtex_splicing()
MIN_STDEV = 0.01


def gtex_splicing(raw_path):
    """

    General handler for all GTEx splicing files.

    Args:
        raw_path (str): the complete path to the
                        raw downloaded file

    Returns:
        Processed DataFrame

    """

    df = pd.read_csv(raw_path, sep="\t")

    df["exon_id"] = concat_cols(
        df,
        [
            "gene_name",
            "event_type",
            "event_chr",
            "event_coordinates",
            "alt_region_coordinates",
        ],
        "_",
    )

    df = df.drop(
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

    df = df.set_index("exon_id")

    gtex_manifest = access.load("gtex_manifest")
    sra_to_gtex = dict(
        zip(gtex_manifest["Comment[ENA_RUN]"], gtex_manifest["Source Name"])
    )

    df.columns = [sra_to_gtex[x[:-8]] for x in df.columns]

    df = df.T
    df = df.astype(np.float16)

    nan_counts = df.isna().sum(axis=0)

    keep_cols = df.columns[nan_counts < len(df) - MIN_VALID_COUNT]

    df = df[keep_cols]

    stdevs = df.std(axis=0)
    df = df[df.columns[stdevs >= MIN_STDEV]]

    df = df.astype(np.float16)

    return df


class Processors:
    """

    GTEx dataset processors.

    """

    @staticmethod
    def gtex_2919_manifest(raw_path):
        """

        Process the GTEx 2919 manifest.

        Args:
            raw_path (str): the complete path to the
                            raw downloaded file

        Returns:
            Processed DataFrame

        """

        df = pd.read_csv(raw_path, sep="\t")

        return df

    @staticmethod
    def gtex_5214_manifest(raw_path):
        """

        Process the GTEx 5214 manifest.

        Args:
            raw_path (str): the complete path to the
                            raw downloaded file

        Returns:
            Processed DataFrame

        """

        df = pd.read_csv(raw_path, sep="\t")

        return df

    @staticmethod
    def gtex_manifest():
        """

        Process the merged manifest.

        Returns:
            Processed DataFrame

        """

        gtex_manifest_1 = access.load("gtex_2919_manifest")
        gtex_manifest_2 = access.load("gtex_5214_manifest")

        gtex_manifest = pd.concat([gtex_manifest_1, gtex_manifest_2], axis=0, sort=True)

        gtex_manifest = gtex_manifest.astype(str)

        return gtex_manifest

    @staticmethod
    def gtex_gene_tpm(raw_path):
        """

        Process the merged manifest.

        Returns:
            Processed DataFrame

        """

        df = pd.read_csv(raw_path, skiprows=2, index_col=0, sep="\t")

        df.index = df["Description"] + "_" + df.index
        df.drop(["Description"], axis=1, inplace=True)
        df = df.T
        df = np.log2(df + 1)

        df = df.astype(np.float16)

        return df

    @staticmethod
    def gtex_a3ss(raw_path):
        """

        Process GTEx A3SS splicing events.

        Args:
            raw_path (str): the complete path to the
                            raw downloaded file

        Returns:
            Processed DataFrame

        """

        return gtex_splicing(raw_path)

    @staticmethod
    def gtex_a5ss(raw_path):
        """

        Process GTEx A5SS splicing events.

        Args:
            raw_path (str): the complete path to the
                            raw downloaded file

        Returns:
            Processed DataFrame

        """

        return gtex_splicing(raw_path)

    @staticmethod
    def gtex_se(raw_path):
        """

        Process GTEx SE splicing events.

        Args:
            raw_path (str): the complete path to the
                            raw downloaded file

        Returns:
            Processed DataFrame

        """

        return gtex_splicing(raw_path)

    @staticmethod
    def gtex_ir(raw_path):
        """

        Process GTEx IR splicing events.

        Args:
            raw_path (str): the complete path to the
                            raw downloaded file

        Returns:
            Processed DataFrame

        """

        return gtex_splicing(raw_path)

    @staticmethod
    def gtex_mx(raw_path):
        """

        Process GTEx MX splicing events.

        Args:
            raw_path (str): the complete path to the
                            raw downloaded file

        Returns:
            Processed DataFrame

        """

        return gtex_splicing(raw_path)
