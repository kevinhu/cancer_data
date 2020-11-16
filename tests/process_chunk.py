import numpy as np
import pandas as pd
from tqdm import tqdm


def row_count(file):
    with open(file) as f:
        for i, l in enumerate(f):
            pass
    return i


def chunk_to_mmap(input_csv):

    columns = pd.read_csv(input_csv, index_col=0, nrows=0).columns
    col_dtypes = {x: np.float16 for x in columns}

    chunk_size = 1000

    chunk_iterator = pd.read_csv(
        input_csv,
        chunksize=chunk_size,
        dtype=col_dtypes,
        engine="c",
        low_memory=False,
        index_col=0,
    )

    n_rows = row_count(input_csv)
    n_cols = len(columns)

    mmap = np.memmap(
        f"{input_csv}.mm",
        dtype="float16",
        mode="w+",
        shape=(n_rows, n_cols),
    )

    chunk_idx = 0

    for chunk in tqdm(chunk_iterator, total=n_rows // chunk_size + 1):

        first_row = chunk_idx * chunk_size
        last_row = (chunk_idx + 1) * chunk_size

        mmap[first_row:last_row] = chunk.values

        chunk_idx += 1

    mmap.flush()


def get_csv_index(file):
    with open(file) as f:
        index = []
        for i, l in enumerate(f):
            index.append(l.partition(",")[0])

        return index[1:]


def mmap_to_hdf(input_csv):

    columns = pd.read_csv(input_csv, index_col=0, nrows=0).columns
    index = get_csv_index(input_csv)

    mmap = np.memmap(
        f"{input_csv}.mm",
        dtype="float16",
        mode="r",
        shape=(len(index), len(columns)),
    )

    df = pd.DataFrame(data=mmap, index=index, columns=columns, copy=False)
    df = df.T
    df.index = df.index.map(lambda x: x.split(".")[0])

    df.to_hdf(
        f"{input_csv}.h5",
        key=str(input_csv),
        complevel=9,
        complib="bzip2",
        mode="w",
    )
