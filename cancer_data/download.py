from .checks import is_downloadable
from .config import DOWNLOAD_DIR, REFERENCE_DIR, SCHEMA
from .utils import download_from_url


def download(dataset_id):
    """

    Handler for downloading a dataset.

    Args:
        dataset_id (str): ID of the dataset

    """

    assert dataset_id in SCHEMA.index, f"{dataset_id} is not in the schema."

    dataset_row = SCHEMA.loc[dataset_id]

    dataset_type = dataset_row["type"]
    dataset_url = dataset_row["url"]
    dataset_downloaded_name = dataset_row["downloaded_name"]
    dataset_downloaded_md5 = dataset_row["downloaded_md5"]

    if dataset_type == "reference":

        download_from_url(
            dataset_url,
            f"{REFERENCE_DIR}/{dataset_downloaded_name}",
            reference_md5=dataset_downloaded_md5,
        )

    elif dataset_type == "primary_dataset":
        download_from_url(
            dataset_url,
            f"{DOWNLOAD_DIR}/{dataset_downloaded_name}",
            reference_md5=dataset_downloaded_md5,
        )

    else:

        raise ValueError("Unsupported dataset type.")


def download_all():
    """

    Download all datasets in the schema.

    """

    for _, dataset in SCHEMA.iterrows():

        dataset_id = dataset["id"]

        if is_downloadable(dataset_id):

            download(dataset_id)
