import requests
import os
import tqdm
from tqdm import tqdm


def download_from_url(url, output_path, overwrite=False):

    if os.path.isfile(output_path):

        if not overwrite:

            print(f"{output_path} already exists, skipping")

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

    return True
