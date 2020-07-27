from .config import DOWNLOAD_DIR, PROCESSED_DIR, PREVIEW_DIR, SCHEMA
from .access import Datasets
from .utils import bcolors, file_exists, export_hdf

from .processors import ccle, depmap, gtex, other, tcga


class Processors(
    ccle.Processors,
    depmap.Processors,
    gtex.Processors,
    other.Processors,
    tcga.Processors,
):
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

    if dependencies is None or dependencies != dependencies:
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

    df = Datasets.load(dataset_id, stop=PREVIEW_LEN)

    df.to_csv(f"{PREVIEW_DIR}/{dataset_id}.txt", sep="\t")


def process(dataset_id, overwrite=False):
    """

    Handler for processing a dataset

    Args:
        dataset_id (str): ID of the dataset
        downloaded_name (str): name of the raw downloaded dataset file
        dependencies (str): dependencies of the dataset
        dataset_type (str): type of the dataset

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

            else:

                print(
                    f"Handler for {id_bold} {bcolors.FAIL}not found{bcolors.ENDC}, skipping"
                )


def process_all():
    """

    Process all datasets in the schema.

    """

    for _, file in SCHEMA.iterrows():

        process(file["id"])
