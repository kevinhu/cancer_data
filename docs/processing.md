## Processing

`cancer_data.check_dependencies(dataset_id)`: Checks if all the dependencies of a dataset are present.

`cancer_data.process(dataset_id, overwrite=False, delete_raw=False)`: Process a downloaded dataset.

`cancer_data.process_all(dataset_id, process_kwargs={})`: Process all datasets.

`cancer_data.download_and_process(dataset_id, download_kwargs={}, process_kwargs={})`: Download and process a dataset.

`cancer_data.download_and_process_all(dataset_id, download_kwargs={}, process_kwargs={})`: Download and process all datasets.
