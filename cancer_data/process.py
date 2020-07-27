import os
import warnings

from . import access
from .config import DOWNLOAD_DIR, PROCESSED_DIR, PREVIEW_DIR, SCHEMA
from .utils import bcolors, file_exists, export_hdf

from .download import download

from .processors import ccle, depmap, gtex, other, tcga


class Processors(
    ccle.Processors,
    depmap.Processors,
    gtex.Processors,
    other.Processors,
    tcga.Processors,
):
    """

    Subclass for merging the processing methods from
    the individual collections.

    """

    def __init__(self):
        return


# row to generate in generate_preview()
PREVIEW_LEN = 10


def check_dependencies(dependencies):
    """

    Check if dataset dependencies are all met.

    Args:
        dependencies (str or NaN): comma-delimited depenencies

    """

    if dependencies is None or dependencies != dependencies or dependencies == "":
        return

    for d in dependencies.split(","):

        d_file = f"{PROCESSED_DIR}/{d}.h5"

        assert file_exists(d_file), f"Dependency {d} does not exist."


def generate_preview(dataset_id):
    """

    Generate a preview of a DataFrame, saving
    to a CSV

    Args:
        dataset_id: the ID of the dataset

    """

    df = access.load(dataset_id, stop=PREVIEW_LEN)

    df.to_csv(f"{PREVIEW_DIR}/{dataset_id}.txt", sep="\t")


def remove_raw(dataset_id):
    """

    Remove the raw dataset file.

    Args:
        dataset_id: ID of dataset to remove

    """

    id_bold = f"{bcolors.BOLD}{dataset_id}{bcolors.ENDC}"

    assert dataset_id in SCHEMA.index, f"{id_bold} is not in the schema."

    dataset_row = SCHEMA.loc[dataset_id]
    downloaded_name = dataset_row["downloaded_name"]

    raw_file = f"{DOWNLOAD_DIR}/{downloaded_name}"

    if file_exists(raw_file):
        os.remove(raw_file)
    else:
        warnings.warn(f"{id_bold} is in schema, but raw file does not exist.")


def remove_processed(dataset_id):
    """

    Remove the processed dataset file.

    Args:
        dataset_id: ID of dataset to remove

    """

    id_bold = f"{bcolors.BOLD}{dataset_id}{bcolors.ENDC}"

    assert dataset_id in SCHEMA.index, f"{id_bold} is not in the schema."

    processed_file = f"{PROCESSED_DIR}/{dataset_id}.h5"

    if file_exists(processed_file):
        os.remove(processed_file)

        print("")
    else:
        warnings.warn(f"{id_bold} is in schema, but processed file does not exist.")


def remove(dataset_id):
    """

    Remove the raw and processed dataset files.

    Args:
        dataset_id: ID of dataset to remove

    """

    assert dataset_id in SCHEMA.index, f"{dataset_id} is not in the schema."

    remove_raw(dataset_id)
    remove_processed(dataset_id)


def remove_all_raw():
    """

    Removes all raw dataset files.

    """

    for _, dataset in SCHEMA.iterrows():

        remove_raw(dataset["id"])


def remove_all_processed():
    """

    Removes all processed dataset files.

    """

    for _, dataset in SCHEMA.iterrows():

        remove_processed(dataset["id"])


def remove_all():
    """

    Removes all raw and processed dataset files.

    """

    for _, dataset in SCHEMA.iterrows():

        remove(dataset["id"])


def process(dataset_id, overwrite=False, delete_raw=False):
    """

    Handler for processing a dataset.

    Args:
        dataset_id (str): ID of the dataset
        overwrite (bool): overwrite existing
        remove_raw (bool): remove the raw file

    """

    assert dataset_id in SCHEMA.index, f"{dataset_id} is not in the schema."

    dataset_row = SCHEMA.loc[dataset_id]

    downloaded_name = dataset_row["downloaded_name"]
    dependencies = dataset_row["dependencies"]
    dataset_type = dataset_row["type"]

    if dataset_type in ["primary_dataset", "secondary_dataset"]:

        output_path = f"{PROCESSED_DIR}/{dataset_id}.h5"

        id_bold = f"{bcolors.BOLD}{dataset_id}{bcolors.ENDC}"

        if file_exists(output_path) and not overwrite:

            print(f"{id_bold} already processed, skipping")

            return

        else:

            handler = getattr(Processors, dataset_id, None)

            if handler is not None:

                print(f"Processing {id_bold}")

                check_dependencies(dependencies)

                if dataset_type in ["primary_dataset"]:
                    df = handler(DOWNLOAD_DIR / downloaded_name)

                elif dataset_type in ["secondary_dataset"]:
                    df = handler()

                export_hdf(dataset_id, df)

                generate_preview(dataset_id)

                if delete_raw:

                    remove_raw(dataset_id)

            else:

                print(
                    f"Handler for {id_bold} {bcolors.FAIL}not found{bcolors.ENDC}, skipping"
                )

                return


def download_and_process(dataset_id, download_kwargs={}, process_kwargs={}):
    """

    Download and process a dataset.

    Args:
        dataset_id (str): ID of the dataset
        download_kwargs (dict): arguments to pass to download()
        process_kwargs (dict): arguments to pass to process()

    """

    id_bold = f"{bcolors.BOLD}{dataset_id}{bcolors.ENDC}"

    assert dataset_id in SCHEMA.index, f"{id_bold} is not in the schema."

    download(dataset_id, **download_kwargs)
    process(dataset_id, **process_kwargs)


def download_and_process_all(dataset_id, download_kwargs={}, process_kwargs={}):
    """

    Download and process all datasets in the schema.

    """

    for _, dataset in SCHEMA.iterrows():

        download_and_process(dataset["id"], download_kwargs, process_kwargs)


def process_all():
    """

    Process all datasets in the schema.

    """

    for _, dataset in SCHEMA.iterrows():

        process(dataset["id"])
