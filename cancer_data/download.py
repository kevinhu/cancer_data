from config import DOWNLOAD_DIR, REFERENCE_DIR, SCHEMA

from utils import download_from_url

if __name__ == "__main__":

    for _, file in SCHEMA.iterrows():

        if file["type"] == "reference":

            download_from_url(
                file["url"],
                f"{REFERENCE_DIR}/{file['downloaded_name']}",
                reference_md5 = file['downloaded_md5']
            )

        elif file["type"] == "primary_dataset":
            download_from_url(
                file["url"],
                f"{DOWNLOAD_DIR}/{file['downloaded_name']}",
                reference_md5 = file['downloaded_md5']
            )