import os
import pathlib
import subprocess
from glob import glob
from typing import Dict, List, Optional

from ffm.acquisition.teleseismicquery import TeleseismicQuery


class CwbQuery:
    def __init__(
        self,
        eventid: str,
        source: str,
    ):
        self.query = TeleseismicQuery.from_id(eventid, source)

    def get_data(
        self,
        directory: pathlib.Path,
        edge_cwb_jar_path: pathlib.Path,
        java: pathlib.Path,
        host: str,
        networks: List[str] = [],
        stations: Optional[Dict[str, Dict[str, List[str]]]] = None,
    ):
        cwb_cmd = f"{str(java)} -jar {str(edge_cwb_jar_path)}"
        commands: List[str] = []
        for network in networks:
            commands += [
                (
                    f"{cwb_cmd} -s {network:.<7}BH. -h {host} "
                    f"-delazc {self.query.min_distance}:{self.query.max_distance}:{self.query.latitude}:{self.query.longitude} "
                    f'-b "{self.date_string}" -d {self.query.duration} -nogaps -sacpz nm'
                )
            ]
        if stations is not None:
            for network_key, network_stations in stations.items():
                for station_key, channels in network_stations.items():
                    for channel in channels:
                        commands += [
                            (
                                f"{cwb_cmd} -s {network_key:.<2}{station_key:.<5}{channel:.<3} -h {host} "
                                f"-delazc {self.query.min_distance}:{self.query.max_distance}:{self.query.latitude}:{self.query.longitude} "
                                f'-b "{self.date_string}" -d {self.query.duration} -nogaps -sacpz nm'
                            )
                        ]
        current_dir = pathlib.Path()
        os.chdir(directory)
        try:
            for command in commands:
                print(f"Querying EDGE CWB: {command}")
                process = subprocess.Popen(command, shell=True)

                stdout, stderr = process.communicate()
                if stdout is not None:
                    print(stdout.decode())
                if process.returncode != 0:
                    if stderr is not None:
                        print("Command error!\n", stderr.decode())
                    else:
                        print("No data")
        finally:
            os.chdir(current_dir)

        # reformat sac PZ files
        for pz in glob("./*.sac.pz"):
            with open(pz, "r") as old:
                old_pz_data = old.read()
            if "* INPUT UNIT   NM" not in old_pz_data:
                print(f"Unit in {pz} not the expected NM. Skipping reformatting")
                continue
            # remove comment lines
            old_lines = old_pz_data.split("\n")
            # remove comments and constant lines
            new_lines = []
            constant: Optional[str] = None
            station: Optional[str] = None
            network_line: Optional[str] = None
            component: Optional[str] = None
            location: Optional[str] = None
            for line in old_lines:
                if line.startswith("*"):
                    if line.startswith("* NETWORK"):
                        network_line = line.split("NETWORK")[-1].strip()
                    elif line.startswith("* STATION"):
                        station = line.split("STATION")[-1].strip()
                    if line.startswith("* COMPONENT"):
                        component = line.split("COMPONENT")[-1].strip()
                    if line.startswith("* LOCATION"):
                        location = line.split("LOCATION")[-1].strip()
                    continue
                elif line.startswith("CONSTANT"):
                    constant = line
                else:
                    new_lines += [line]
            if (
                constant is None
                or station is None
                or network_line is None
                or component is None
                or location is None
            ):
                print(
                    f"CONSTANT and/or station information cannot be found in {pz}. Skipping reformatting"
                )
                continue
            constant_val = constant.split("CONSTANT")[-1].strip()
            new_constant = float(constant_val) * 1e9
            updated_constant = constant.replace(constant_val, f"{new_constant:.4e}")
            if new_lines[-1] == "":
                new_lines[-1] = updated_constant
            else:
                new_lines += [updated_constant]
            new_file = f"SAC_PZs_{network_line}_{station}_{component}_{location}"
            with open(new_file, "w") as new:
                new.write("\n".join(new_lines))
            os.remove(pz)
        return commands

    @property
    def date_string(self) -> str:
        """Format date string for CWB query"""
        return self.query.event_time.strftime("%Y/%m/%d %H:%M:%S")
