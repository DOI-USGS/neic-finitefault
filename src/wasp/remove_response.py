import glob
import os
import pathlib
from typing import List, Union

from obspy.io.sac import SACTrace  # type: ignore


def get_sacpz_file(sacheader: SACTrace, data: str = "obspy") -> str:
    """Name the sacpz file

    :param sacheader: The sac trace header
    :type sacheader: SACTrace
    :param data: The data type, defaults to "obspy"
    :type data: str, optional
    :return: The file name
    :rtype: str
    """
    loc = sacheader.khole
    if data == "IRIS" and not loc:
        loc = "--"
    else:
        if not loc:
            loc = "__"
        elif len(loc) == 0:
            loc = "__"
    pzfile = "SAC_PZs_{}_{}_{}_{}".format(
        sacheader.knetwk, sacheader.kstnm, sacheader.kcmpnm, loc
    )
    if data == "IRIS":
        pzfile = "SACPZ.{}.{}.{}.{}".format(
            sacheader.knetwk, sacheader.kstnm, loc, sacheader.kcmpnm
        )
    return pzfile


def convert_response_acc(
    resp_file: Union[pathlib.Path, str]
) -> Union[pathlib.Path, str]:
    """Remove zeros from the response file

    :param resp_file: The response file location
    :type resp_file: Union[pathlib.Path, str]
    :return: The response file location
    :rtype: Union[pathlib.Path, str]
    """
    with open(resp_file, "r") as infile:
        lines = [line for line in infile]
    indexes = [i for i, line in enumerate(lines) if "ZEROS" in line.split()]
    index0 = 0
    lines2: List[str] = []
    for index in indexes:
        lines2 = lines2 + lines[index0:index]
        string, nzeroes = lines[index].split()
        if nzeroes == "2":
            lines2 = lines2 + ["ZEROS\t 0\n"]
            index0 = index + 3
        if nzeroes == "3":
            lines2 = lines2 + ["ZEROS\t 0\n"]
            index0 = index + 4
        # TODO: Investigate: What is happening here. Is it just removing the zeros?
        # What about ones with zeros greater than 3?
        lines2 = lines2 + lines[index0:]
    with open(resp_file, "w") as outfile:
        for line in lines2:
            outfile.write(line)
    return resp_file


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f", "--folder", default=os.getcwd(), help="folder where there are input files"
    )
    args = parser.parse_args()
    os.chdir(args.folder)
    tensor_info = {"depth": 25, "time_shift": 40}
    data_prop = {
        "tele_filter": {"freq0": 0.001, "freq1": 0.002, "freq2": 1.0, "freq3": 1.2}
    }
    files0 = glob.glob("SACPZ*")
    for file in files0:
        convert_response_acc(file)
    files = glob.glob("*SAC")
    # __remove_response_str2(files, data="IRIS")
