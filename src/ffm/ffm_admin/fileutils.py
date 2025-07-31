import pathlib
from typing import List


def validate_files(files: List[pathlib.Path]):
    "Validate files exist"
    for file in files:
        file = pathlib.Path(file)
        if not file.exists():
            raise FileNotFoundError(
                f"Expected file ({str(file.resolve())}) does not exist!"
            )
