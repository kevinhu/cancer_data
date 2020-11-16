import os
import pathlib

import pandas as pd

"""
Contains configurations and settings used by the rest of the project.
"""

_PROJECT_DIR = pathlib.Path(__file__).resolve().parent

DATA_DIR = _PROJECT_DIR / "data"
DOWNLOAD_DIR = DATA_DIR / "raw"
REFERENCE_DIR = DATA_DIR / "reference"
PROCESSED_DIR = DATA_DIR / "processed"
PREVIEW_DIR = DATA_DIR / "preview"

SCHEMA = pd.read_csv(_PROJECT_DIR / "schema.csv")
SCHEMA.index = SCHEMA["id"]
SCHEMA["dependencies"] = SCHEMA["dependencies"].fillna("")
