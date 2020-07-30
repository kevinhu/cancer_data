from .config import SCHEMA


def is_downloadable(dataset_id):
    """

    Check if a dataset is downloadable.

    Args:
        dataset_id (str): ID of the dataset

    """

    dataset_row = SCHEMA.loc[dataset_id]
    dataset_type = dataset_row["type"]

    if dataset_type in ["reference", "primary_dataset"]:

        return True

    return False


def is_processable(dataset_id):
    """

    Check if a dataset can be processed.

    Args:
        dataset_id (str): ID of the dataset

    """

    dataset_row = SCHEMA.loc[dataset_id]
    dataset_type = dataset_row["type"]

    if dataset_type in ["primary_dataset", "secondary_dataset"]:

        return True

    return False
