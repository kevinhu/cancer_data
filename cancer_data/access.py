import pandas as pd

from .checks import is_downloadable, is_processable
from .config import DOWNLOAD_DIR, PREVIEW_DIR, PROCESSED_DIR, REFERENCE_DIR, SCHEMA
from .utils import bcolors, file_exists


def load(dataset_id, **read_hdf_kwargs):
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

    df = pd.read_hdf(dataset_path, **read_hdf_kwargs)

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


def status():
    """

    Print the statuses of all datasets in the schema.

    """

    max_id_len = SCHEMA["id"].apply(len).max()

    for _, dataset in SCHEMA.iterrows():

        dataset_id = dataset["id"]
        dataset_type = dataset["type"]
        downloaded_name = dataset["downloaded_name"]

        id_bold = f"{bcolors.BOLD}{dataset_id}{bcolors.ENDC}"

        print(f"{id_bold}:", end="")
        print((max_id_len - len(dataset_id)) * " ", end="")

        if is_downloadable(dataset_id):

            download_path = DOWNLOAD_DIR / downloaded_name

            if dataset_type == "reference":

                download_path = REFERENCE_DIR / downloaded_name

            if file_exists(download_path):
                print(f"{bcolors.OKGREEN}downloaded{bcolors.ENDC}\t", end="")
            else:
                print(f"{bcolors.FAIL}not downloaded{bcolors.ENDC}\t", end="")

        else:

            print("\t\t", end="")

        if is_processable(dataset_id):
            if file_exists(PROCESSED_DIR / f"{dataset_id}.h5"):
                print(f"{bcolors.OKGREEN}processed{bcolors.ENDC}", end="")
            else:
                print(f"{bcolors.FAIL}not processed{bcolors.ENDC}", end="")

        print()


def summary(dataset_id):
    """

    Print the summary info of a dataset.

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

    print(dataset_summary)


def schema():
    """

    Return a copy of the schema

    """

    return SCHEMA.copy(deep=True)


def types():
    """

    Return a sorted list of all the dataset types

    """

    return sorted(list(set(SCHEMA["type"])))


def sources():
    """

    Return a sorted list of all the dataset sources

    """

    print(SCHEMA["source"])

    return sorted(list(set(SCHEMA["source"])))
