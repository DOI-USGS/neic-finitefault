import pathlib
from typing import List


def validate_files(files: List[pathlib.Path]):
    "Validate files exist"
    print("HERE", files)
    for file in files:
        print("The file", file)
        if not file.exists():
            raise FileNotFoundError(
                f"Expected file ({str(file.resolve())}) does not exist!"
            )
