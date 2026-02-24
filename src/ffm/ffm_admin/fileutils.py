import pathlib
from typing import List, Optional


def validate_files(files: List[pathlib.Path]):
    "Validate files exist"
    for file in files:
        file = pathlib.Path(file)
        if not file.exists():
            raise FileNotFoundError(
                f"Expected file ({str(file.resolve())}) does not exist!"
            )


def validate_data_directory(directory: Optional[pathlib.Path]) -> pathlib.Path:
    """Validate directory or create data directory"""
    if directory is None:
        directory = pathlib.Path() / "data"
        directory.mkdir(exist_ok=True)

    if not directory.exists():
        raise FileNotFoundError(f"Directory ({directory}) does not exist")

    return directory
