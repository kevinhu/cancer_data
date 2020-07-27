from cancer_data import process, download, access

df = access.Datasets.load("ccle_annotations")

print(df)

