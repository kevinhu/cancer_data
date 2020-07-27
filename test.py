import cancer_data

cancer_data.download_and_process("ccle_annotations", process_kwargs={"overwrite": True,"delete_raw":True})

print(cancer_data.summary("ccle_annotations"))
