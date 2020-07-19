import os
import pathlib
import pandas as pd

_PROJECT_DIR = pathlib.Path(__file__).resolve().parent

DATA_DIR = _PROJECT_DIR / "data"
DOWNLOAD_DIR = DATA_DIR / "raw"
REFERENCE_DIR = DATA_DIR / "reference"
PROCESSED_DIR = DATA_DIR / "processed"

SCHEMA = pd.read_csv(_PROJECT_DIR / "schema.txt", sep="\t")
