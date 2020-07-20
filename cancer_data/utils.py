import requests
import os
import tqdm
from tqdm import tqdm
from hashlib import md5


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
    with open(file_path, "rb") as f:

        data = f.read()

        file_md5 = md5(data).hexdigest()

        return file_md5 == reference_md5


def download_from_url(
    url, output_path, overwrite=False, reference_md5=None, is_retry=False
):

    if os.path.isfile(output_path):

        if not overwrite:

            if reference_md5 is not None:

                if not md5_match(output_path, reference_md5):

                    if not is_retry:

                        print(
                            f"{output_path} {bcolors.FAIL}does not match{bcolors.ENDC} provided md5sum. Attempting download."
                        )

                else:

                    print(
                        f"{output_path} already downloaded and {bcolors.OKGREEN}matches{bcolors.ENDC} provided md5sum."
                    )

                    return True

            else:

                print(
                    f"{output_path} already exists, skipping md5sum check (not provided)"
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
                    f"{output_path} {bcolors.FAIL}does not match{bcolors.ENDC} provided md5sum. Attempting second download."
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
                    f"{bcolors.FAIL}Second download of {output_path} failed. Recommend manual inspection.{bcolors.ENDC}"
                )

                return False

        else:

            print(
                f"{output_path} {bcolors.OKGREEN}matches{bcolors.ENDC} provided md5sum."
            )

    return True
