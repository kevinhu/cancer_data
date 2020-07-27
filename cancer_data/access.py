import pandas as pd

from .config import DOWNLOAD_DIR, PROCESSED_DIR, PREVIEW_DIR, SCHEMA
from .utils import file_exists


class Datasets:
    """

    General class for interacting with datasets.

    """

    @staticmethod
    def load(dataset_id, **kwargs):
        """

        Load a processed dataset.

        Args:
            dataset_id (str): ID of the dataset
            **kwards: additional arguments to pass to pd.read_hdf()

        """

        dataset_path = PROCESSED_DIR / f"{dataset_id}.h5"

        assert dataset_id in SCHEMA.index, f"{dataset_id} not in schema!"
        assert file_exists(dataset_path), f"{dataset_id} does not exist!"

        df = pd.read_hdf(dataset_path, **kwargs)

        return df

    @staticmethod
    def description(dataset_id):
        """

        Get the description of a dataset.

        Args:
            dataset_id (str): ID of the dataset

        """

        assert dataset_id in SCHEMA.index, f"{dataset_id} not in schema!"

        return SCHEMA.loc[dataset_id, "description"]

    @staticmethod
    def summary(dataset_id):
        """

        Get the summary info of a dataset.

        Args:
            dataset_id (str): ID of the dataset

        """

        assert dataset_id in SCHEMA.index, f"{dataset_id} not in schema!"

        dataset_row = SCHEMA.loc[dataset_id]

        if dataset_row["type"] == "reference":

            dataset_summary = "\n".join(
                [
                    f"ID: {dataset_row['id']}",
                    f"Type: {dataset_row['type']}",
                    f"Description: {dataset_row['description']}",
                    f"Provider: {dataset_row['source']}",
                    f"Source URL: {dataset_row['url']}",
                    f"Portal URL: {dataset_row['portal_url']}",
                    f"Size: {dataset_row['downloaded_size']}",
                    f"md5sum: {dataset_row['downloaded_md5']}",
                ]
            )

        elif dataset_row["type"] == "primary_dataset":

            dataset_summary = "\n".join(
                [
                    f"ID: {dataset_row['id']}",
                    f"Type: {dataset_row['type']}",
                    f"Dependencies: {dataset_row['dependencies']}",
                    f"Description: {dataset_row['description']}",
                    f"Provider: {dataset_row['source']}",
                    f"Source URL: {dataset_row['url']}",
                    f"Portal URL: {dataset_row['portal_url']}",
                    f"Size: {dataset_row['downloaded_size']}",
                    f"md5sum: {dataset_row['downloaded_md5']}",
                ]
            )

        elif dataset_row["type"] == "secondary_dataset":

            dataset_summary = "\n".join(
                [
                    f"ID: {dataset_row['id']}",
                    f"Type: {dataset_row['type']}",
                    f"Dependencies: {dataset_row['dependencies']}",
                    f"Description: {dataset_row['description']}",
                    f"Provider: {dataset_row['source']}",
                    f"Source URL: {dataset_row['url']}",
                    f"Portal URL: {dataset_row['portal_url']}",
                    f"Size: {dataset_row['downloaded_size']}",
                    f"md5sum: {dataset_row['downloaded_md5']}",
                ]
            )

        return dataset_summary
