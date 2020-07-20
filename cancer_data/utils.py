import requests
import os
import tqdm
from tqdm import tqdm
from hashlib import md5

from pathlib import Path


class bcolors:
    HEADER = "\033[95m"
    OKBLUE = "\033[94m"
    OKGREEN = "\033[92m"
    WARNING = "\033[93m"
    FAIL = "\033[91m"
    ENDC = "\033[0m"
    BOLD = "\033[1m"
    UNDERLINE = "\033[4m"


def md5_match(file_path, reference_md5):

    """
    Checks if a file matches a provided md5sum

    Parameters
    ----------
    file_path: string
        path to the file to check
    reference_md5: string
        md5sum to check

    Returns
    -------
    True if the file matches the md5sum, False otherwise

    """

    with open(file_path, "rb") as f:

        data = f.read()

        file_md5 = md5(data).hexdigest()

        return file_md5 == reference_md5


def download_from_url(
    url, output_path, overwrite=False, reference_md5=None, is_retry=False
):

    """
    Download a file from a URL.

    Shows progress bar and checks md5sum. Also 
    checks if file is already downloaded.

    Parameters
    ----------
    url: string
        URL to download from
    output_path: string
        path to save the output to
    overwrite: boolean
        whether or not to overwrite the file if it already exists
    reference_md5: string
        md5sum to check
    is_retry:
        whether or not the download is a retry (if the md5sum)
        does not match, is called again

    Returns
    -------
    True if download was successful, False otherwise

    """

    output_filename = Path(output_path).name
    output_filename_bold = f"{bcolors.BOLD}{output_filename}{bcolors.ENDC}"

    if os.path.isfile(output_path):

        if not overwrite:

            if reference_md5 is not None:

                if not md5_match(output_path, reference_md5):

                    if not is_retry:

                        print(
                            f"{output_filename_bold} {bcolors.FAIL}does not match{bcolors.ENDC} provided md5sum. Attempting download."
                        )

                else:

                    print(
                        f"{output_filename_bold} already downloaded and {bcolors.OKGREEN}matches{bcolors.ENDC} provided md5sum."
                    )

                    return True

            else:

                print(
                    f"{output_filename_bold} already exists, skipping md5sum check (not provided)"
                )

                return True

    print(f"Downloading {url}")

    r = requests.get(url, stream=True)

    total_size = int(r.headers.get("content-length", 0))
    block_size = 1024

    t = tqdm(total=total_size, unit="iB", unit_scale=True)

    with open(output_path, "wb") as f:
        for data in r.iter_content(block_size):
            t.update(len(data))
            f.write(data)

    t.close()

    if total_size != 0 and t.n != total_size:
        print("Download error: sizes do not match")

        return False

    if reference_md5 is not None:

        if not md5_match(output_path, reference_md5):

            if not is_retry:

                print(
                    f"{output_filename_bold} {bcolors.FAIL}does not match{bcolors.ENDC} provided md5sum. Attempting second download."
                )

                download_from_url(
                    url,
                    output_path,
                    overwrite=overwrite,
                    reference_md5=reference_md5,
                    is_retry=True,
                )

            else:

                print(
                    f"{bcolors.FAIL}Second download of {output_filename_bold} failed. Recommend manual inspection.{bcolors.ENDC}"
                )

                return False

        else:

            print(
                f"{output_filename_bold} {bcolors.OKGREEN}matches{bcolors.ENDC} provided md5sum."
            )

    return True
