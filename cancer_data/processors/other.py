from gtfparse import read_gtf


class Processors:
    """

    Other dataset processors.

    """

    @staticmethod
    def g19_7_definitions(raw_path):
        """

        Process GENCODE g19 v7 gene+transcript definitions.

        Args:
            raw_path (str): the complete path to the
                            raw downloaded file

        Returns:
            Processed DataFrame

        """

        df = read_gtf(raw_path)

        return df

    @staticmethod
    def ensembl_75_definitions(raw_path):
        """

        Process ENSEMBL v75 gene_transcript definitions.

        Args:
            raw_path (str): the complete path to the
                            raw downloaded file

        Returns:
            Processed DataFrame

        """

        df = read_gtf(raw_path)

        return df
