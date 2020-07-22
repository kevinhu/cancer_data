from config import DOWNLOAD_DIR, PROCESSED_DIR, SCHEMA

import os

import numpy as np
import pandas as pd

from utils import bcolors, file_exists, export_hdf

from gtfparse import read_gtf

from collections import defaultdict

def check_dependencies(dependencies):

    for d in dependencies.split(","):

        d_file = f"{PROCESSED_DIR}/{d}.h5"

        assert file_exists(d_file), f"Dependency {d} does not exist."


class Processors:
    def __init__(self):
        return

    def g19_7_definitions(raw_path, output_id, dependencies=None):
        df = read_gtf(raw_path)

        export_hdf(output_id, df)

    def ensembl_75_definitions(raw_path, output_id, dependencies=None):
        df = read_gtf(raw_path)

        export_hdf(output_id, df)

    def gtex_2919_manifest(raw_path, output_id, dependencies=None):
        df = pd.read_csv(raw_path, sep="\t")

        export_hdf(output_id, df)

    def gtex_5214_manifest(raw_path, output_id, dependencies=None):
        df = pd.read_csv(raw_path, sep="\t")

        export_hdf(output_id, df)

    def gtex_manifest(raw_path, output_id, dependencies):

        check_dependencies(dependencies)

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

    # TODO: GTEx splicing

    def ccle_annotations(raw_path, output_id, dependencies=None):

        df = pd.read_csv(raw_path, sep="\t")
        df = df.astype(str)

        export_hdf(output_id, df)

    def ccle_translocations_svaba(raw_path, output_id, dependencies):

        check_dependencies(dependencies)

        df = pd.read_excel(raw_path)

        ccle_annotations = pd.read_hdf(f"{PROCESSED_DIR}/ccle_annotations.h5")
        ccle_to_depmap = dict(
            zip(ccle_annotations["CCLE_ID"], ccle_annotations["depMapID"])
        )

        string_cols = [
            "CCLE_name",
            "map_id",
            "bp1",
            "bp2",
            "class",
            "gene1",
            "gene2",
            "site1",
            "site2",
            "fusion",
        ]

        for col in string_cols:
            df[col] = df[col].astype(str)

        df["depmap_id"] = df["CCLE_name"].apply(ccle_to_depmap.get)

        export_hdf(output_id, df)

    def ccle_rppa_info(raw_path, output_id, dependencies=None):

        df = pd.read_csv(raw_path)
        df = df.astype(str)

        df["format_id"] = (
            df["Target_Genes"].apply(lambda x: x.replace(" ", "-"))
            + "_"
            + df["Antibody_Name"]
        )

        export_hdf(output_id, df)

    def ccle_rppa(raw_path, output_id, dependencies):

        check_dependencies(dependencies)

        df = pd.read_csv(raw_path, index_col=0)

        ccle_rppa_info = pd.read_hdf(f"{PROCESSED_DIR}/ccle_rppa_info.h5")
        antibody_name_map = dict(
            zip(ccle_rppa_info["Antibody_Name"], ccle_rppa_info["format_id"])
        )

        ccle_annotations = pd.read_hdf(f"{PROCESSED_DIR}/ccle_annotations.h5")
        ccle_to_depmap = dict(
            zip(ccle_annotations["CCLE_ID"], ccle_annotations["depMapID"])
        )

        df.columns = map(antibody_name_map.get, df.columns)
        df.index = df.index.map(ccle_to_depmap.get)

        export_hdf(output_id, df)

    def ccle_gene_tpm(raw_path, output_id, dependencies):

        check_dependencies(dependencies)

        df = pd.read_csv(raw_path, sep="\t", index_col=0)
        df = df.iloc[:, 1:]

        g19_7_definitions = pd.read_hdf(f"{PROCESSED_DIR}/g19_7_definitions.h5")

        gene_name_map = dict(zip(g19_7_definitions["gene_id"],g19_7_definitions["gene_name"]))
        gene_name_map = defaultdict(str, gene_name_map)

        gene_names = df.index.map(lambda x: f"{gene_name_map.get(x)}_{x}")

        df = df.set_index(gene_names)
        df = np.log2(df+1)
        df = df.astype(np.float16)
        
        df = df.T

        ccle_annotations = pd.read_hdf(f"{PROCESSED_DIR}/ccle_annotations.h5")
        ccle_to_depmap = dict(
            zip(ccle_annotations["CCLE_ID"], ccle_annotations["depMapID"])
        )

        df.index = df.index.map(ccle_to_depmap.get)

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
