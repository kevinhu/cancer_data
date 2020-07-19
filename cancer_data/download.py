from config import DOWNLOAD_DIR

from utils import download_from_url

# hg19 FASTA 
download_from_url(
    "https://data.broadinstitute.org/snowman/hg19/genomestrip/Homo_sapiens_assembly19.fasta",
    f"{DOWNLOAD_DIR}/Homo_sapiens_assembly19.fasta",
)

# hg19 FASTA index
download_from_url(
    "https://data.broadinstitute.org/snowman/hg19/genomestrip/Homo_sapiens_assembly19.fasta.fai",
    f"{DOWNLOAD_DIR}/Homo_sapiens_assembly19.fasta.fai",
)

# hg19 FASTA dictionary
download_from_url(
    "https://data.broadinstitute.org/snowman/hg19/genomestrip/Homo_sapiens_assembly19.dict",
    f"{DOWNLOAD_DIR}/Homo_sapiens_assembly19.dict",
)
