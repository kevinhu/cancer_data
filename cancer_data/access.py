import pandas as pd

from .config import DOWNLOAD_DIR, PROCESSED_DIR, PREVIEW_DIR, SCHEMA
from .utils import bcolors, file_exists


def load(dataset_id, **kwargs):
    """

    Load a processed dataset.

    Args:
        dataset_id (str): ID of the dataset
        **kwards: additional arguments to pass to pd.read_hdf()

    """

    id_bold = f"{bcolors.BOLD}{dataset_id}{bcolors.ENDC}"

    dataset_path = PROCESSED_DIR / f"{dataset_id}.h5"

    assert dataset_id in SCHEMA.index, f"{id_bold} not in schema."
    assert file_exists(dataset_path), f"{id_bold} does not exist."

    df = pd.read_hdf(dataset_path, **kwargs)

    return df


def description(dataset_id):
    """

    Get the description of a dataset.

    Args:
        dataset_id (str): ID of the dataset

    """

    id_bold = f"{bcolors.BOLD}{dataset_id}{bcolors.ENDC}"

    assert dataset_id in SCHEMA.index, f"{id_bold} not in schema."

    return SCHEMA.loc[dataset_id, "description"]


def summary(dataset_id):
    """

    Get the summary info of a dataset.

    Args:
        dataset_id (str): ID of the dataset

    """

    id_bold = f"{bcolors.BOLD}{dataset_id}{bcolors.ENDC}"

    assert dataset_id in SCHEMA.index, f"{id_bold} not in schema."

    dataset_row = SCHEMA.loc[dataset_id]

    raw_path = DOWNLOAD_DIR / dataset_row["downloaded_name"]
    processed_path = PROCESSED_DIR / f"{dataset_id}.h5"

    yes = f"{bcolors.OKGREEN}yes{bcolors.ENDC}"
    no = f"{bcolors.FAIL}no{bcolors.ENDC}"

    if dataset_row["type"] == "reference":

        dataset_summary = "\n".join(
            [
                f"{bcolors.BOLD}ID{bcolors.ENDC}: {dataset_row['id']}",
                f"{bcolors.BOLD}Type{bcolors.ENDC}: {dataset_row['type']}",
                f"{bcolors.BOLD}Description{bcolors.ENDC}: {dataset_row['description']}",
                f"{bcolors.BOLD}Source{bcolors.ENDC}: {dataset_row['source']}",
                f"{bcolors.BOLD}Source URL{bcolors.ENDC}: {dataset_row['url']}",
                f"{bcolors.BOLD}Portal URL{bcolors.ENDC}: {dataset_row['portal_url']}",
                f"{bcolors.BOLD}Size{bcolors.ENDC}: {int(dataset_row['downloaded_size'])}",
                f"{bcolors.BOLD}md5sum{bcolors.ENDC}: {dataset_row['downloaded_md5']}",
                f"{bcolors.BOLD}Raw path{bcolors.ENDC}: {DOWNLOAD_DIR}/{dataset_row['downloaded_name']}",
                f"{bcolors.BOLD}Raw exists{bcolors.ENDC}: {yes if file_exists(raw_path) else no}",
            ]
        )

    elif dataset_row["type"] == "primary_dataset":

        dataset_summary = "\n".join(
            [
                f"{bcolors.BOLD}ID{bcolors.ENDC}: {dataset_row['id']}",
                f"{bcolors.BOLD}Type{bcolors.ENDC}: {dataset_row['type']}",
                f"{bcolors.BOLD}Dependencies{bcolors.ENDC}: {dataset_row['dependencies']}",
                f"{bcolors.BOLD}Description{bcolors.ENDC}: {dataset_row['description']}",
                f"{bcolors.BOLD}Source{bcolors.ENDC}: {dataset_row['source']}",
                f"{bcolors.BOLD}Source URL{bcolors.ENDC}: {dataset_row['url']}",
                f"{bcolors.BOLD}Portal URL{bcolors.ENDC}: {dataset_row['portal_url']}",
                f"{bcolors.BOLD}Size{bcolors.ENDC}: {int(dataset_row['downloaded_size'])}",
                f"{bcolors.BOLD}md5sum{bcolors.ENDC}: {dataset_row['downloaded_md5']}",
                f"{bcolors.BOLD}Raw path{bcolors.ENDC}: {DOWNLOAD_DIR}/{dataset_row['downloaded_name']}",
                f"{bcolors.BOLD}Raw exists{bcolors.ENDC}: {yes if file_exists(raw_path) else no}",
                f"{bcolors.BOLD}Processed path{bcolors.ENDC}: {PROCESSED_DIR}/{dataset_id}",
                f"{bcolors.BOLD}Processed exists{bcolors.ENDC}: {yes if file_exists(processed_path) else no}",
            ]
        )

    elif dataset_row["type"] == "secondary_dataset":

        dataset_summary = "\n".join(
            [
                f"{bcolors.BOLD}ID{bcolors.ENDC}: {dataset_row['id']}",
                f"{bcolors.BOLD}Type{bcolors.ENDC}: {dataset_row['type']}",
                f"{bcolors.BOLD}Dependencies{bcolors.ENDC}: {dataset_row['dependencies']}",
                f"{bcolors.BOLD}Description{bcolors.ENDC}: {dataset_row['description']}",
                f"{bcolors.BOLD}Source{bcolors.ENDC}: {dataset_row['source']}",
                f"{bcolors.BOLD}Source URL{bcolors.ENDC}: {dataset_row['url']}",
                f"{bcolors.BOLD}Portal URL{bcolors.ENDC}: {dataset_row['portal_url']}",
                f"{bcolors.BOLD}Size{bcolors.ENDC}: {dataset_row['downloaded_size']}",
                f"{bcolors.BOLD}md5sum{bcolors.ENDC}: {dataset_row['downloaded_md5']}",
                f"{bcolors.BOLD}Processed path{bcolors.ENDC}: {PROCESSED_DIR}/{dataset_id}",
                f"{bcolors.BOLD}Processed exists{bcolors.ENDC}: {yes if file_exists(processed_path) else no}",
            ]
        )

    return dataset_summary
