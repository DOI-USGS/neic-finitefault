import pathlib
from copy import deepcopy
from typing import List, Optional


def update_manager_file_locations(
    manage_dicts: List[dict],
    new_directory: pathlib.Path,
    replace_dir: Optional[str] = None,
) -> List[dict]:
    """Update file locations for testing"""
    updated_managed_dicts = deepcopy(manage_dicts)
    for idx, d in enumerate(manage_dicts):
        if replace_dir is not None:
            updated_managed_dicts[idx]["file"] = manage_dicts[idx]["file"].replace(
                replace_dir, str(new_directory)
            )
        else:
            updated_managed_dicts[idx]["file"] = str(
                new_directory / manage_dicts[idx]["file"]
            )
    return updated_managed_dicts
