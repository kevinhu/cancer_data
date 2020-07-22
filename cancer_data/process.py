from config import DOWNLOAD_DIR, PROCESSED_DIR, SCHEMA

import os

import numpy as np
import pandas as pd

from utils import bcolors, file_exists, export_hdf

from gtfparse import read_gtf

class Processors:
    def __init__(self):
        return

    def g19_7_definitions(raw_path, output_id, dependencies=None):
        df = read_gtf(raw_path)

        export_hdf(output_id, df)

    def gtex_2919_manifest(raw_path, output_id, dependencies=None):
        df = pd.read_csv(raw_path, sep="\t")

        export_hdf(output_id, df)

    def gtex_5214_manifest(raw_path, output_id, dependencies=None):
        df = pd.read_csv(raw_path, sep="\t")

        export_hdf(output_id, df)

    def gtex_manifest(raw_path, output_id, dependencies):

        for d in dependencies.split(","):

            d_file = f"{PROCESSED_DIR}/{d}.h5"

            assert file_exists(d_file), f"Dependency {d} does not exist."

        gtex_manifest_1 = pd.read_hdf(f"{PROCESSED_DIR}/gtex_2919_manifest.h5")
        gtex_manifest_2 = pd.read_hdf(f"{PROCESSED_DIR}/gtex_5214_manifest.h5")

        gtex_manifest = pd.concat([gtex_manifest_1, gtex_manifest_2], axis=0, sort=True)

        gtex_manifest = gtex_manifest.astype(str)

        export_hdf(output_id, gtex_manifest)

    def gtex_gene_tpm(raw_path, output_id, dependencies=None):
        df = pd.read_csv(raw_path, skiprows=2, index_col=0, sep="\t")

        df.index = df["Description"] + "_" + df.index
        df.drop(["Description"], axis=1, inplace=True)
        df = df.T
        df = np.log2(df + 1)

        df = df.astype(np.float16)

        export_hdf(output_id, df)




if __name__ == "__main__":

    for _, file in SCHEMA.iterrows():

        if file["type"] in ["primary_dataset", "secondary_dataset"]:

            output_path = f"{PROCESSED_DIR}/{file['id']}.h5"

            if file_exists(output_path):

                print(
                    f"{bcolors.BOLD}{file['id']}{bcolors.ENDC} already processed, skipping"
                )

            else:

                handler = getattr(Processors, file["id"], None)

                if handler is not None:

                    print(f"Processing {bcolors.BOLD}{file['id']}{bcolors.ENDC}")
                    handler(
                        f"{DOWNLOAD_DIR}/{file['downloaded_name']}",
                        file["id"],
                        file["dependencies"],
                    )
