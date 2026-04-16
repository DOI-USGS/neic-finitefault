import os
import pathlib
import subprocess
from typing import Dict, List, Optional

from ffm.acquisition.teleseismicquery import TeleseismicQuery


class CwbQuery(TeleseismicQuery):
    def get_data(
        self,
        directory: pathlib.Path,
        edge_cwb_jar_path: pathlib.Path,
        java: pathlib.Path,
        networks: List[str] = [],
        stations: Optional[Dict[str, Dict[str, List[str]]]] = None,
    ):
        cwb_cmd = f"{str(java)} -jar {str(edge_cwb_jar_path)}"
        commands: List[str] = []
        for network in networks:
            commands += [
                (
                    f"{cwb_cmd} -s {network:.<7}BH. "
                    f"-delazc {self.min_distance}:{self.max_distance}:{self.latitude}:{self.longitude} "
                    f'-b "{self.date_string}" -d {self.duration} -nogaps -sacpz nm'
                )
            ]
        if stations is not None:
            for network_key, network_stations in stations.items():
                for station_key, channels in network_stations.items():
                    for channel in channels:
                        commands += [
                            (
                                f"{cwb_cmd} -s {network_key:.<2}{station_key:.<5}{channel:.<3} "
                                f"-delazc {self.min_distance}:{self.max_distance}:{self.latitude}:{self.longitude} "
                                f'-b "{self.date_string}" -d {self.duration} -nogaps -sacpz nm'
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
        return commands

    @property
    def date_string(self) -> str:
        """Format date string for CWB query"""
        return self.event_time.strftime("%Y/%m/%d %H:%M:%S")
