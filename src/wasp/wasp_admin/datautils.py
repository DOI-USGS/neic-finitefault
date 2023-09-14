from typing import List


def validate_data_types(data_types: List[str], allowed_types: List[str]):
    """Validate data types"""
    for dt in data_types:
        if dt not in allowed_types:
            raise ValueError(
                f"'{dt}' is not in the allowed data type list: {allowed_types}."
            )
