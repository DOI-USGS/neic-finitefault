import json
import logging
import os
import pathlib
from glob import glob
from shutil import move
from typing import List, Tuple

import typer

from wasp.eventpage_downloads import (
    make_waveproperties_json,
    temporary_file_reorganization_for_publishing,
    write_CMTSOLUTION_file,
    write_Coulomb_file,
    write_Okada_displacements,
)
from wasp.fault_plane import point_sources_param, shear_modulous
from wasp.get_outputs import get_insar, read_solution_static_format, retrieve_gps
from wasp.load_ffm_model import load_ffm_model
from wasp.management import default_dirs
from wasp.plot_graphic_NEIC import (
    PlotComparisonMap,
    PlotInsar,
    PlotSlipDist_Compare,
    calculate_cumulative_moment_tensor,
    plot_ffm_sol,
    plot_misfit,
    shakemap_polygon,
)
from wasp.plot_Map import PlotMap
from wasp.read_config import CONFIG_PATH
from wasp.seismic_tensor import get_tensor
from wasp.static2fsp import static_to_fsp
from wasp.static2srf import static_to_srf
from wasp.velocity_models import select_velmodel
from wasp.write_KML import PlotMap_KML

from .datautils import DEFAULT_MANAGEMENT_FILES, ManagedDataTypes
from .fileutils import validate_files

app = typer.Typer(help="Plotting scripts")


@app.command(help="Plot map as KML")
def kml(
    directory: pathlib.Path = typer.Argument(..., help="Path to the data directory"),
    gcmt_tensor_file: pathlib.Path = typer.Argument(
        ..., help="Path to the GCMT moment tensor file"
    ),
    config_file: pathlib.Path = typer.Option(
        CONFIG_PATH, "-c", "--config-file", help="Path to config file"
    ),
    data_types: List[ManagedDataTypes] = typer.Option(
        [],
        "-t",
        "--data-type",
        help="Data types to plot misfit for, default is []",
    ),
    event_id: str = typer.Option(
        None,
        "-e",
        "--event-id",
        help="The event's name/ID",
    ),
    legend_length: float = typer.Option(
        None,
        "-ll",
        "--legend-length",
        help="Set the length (cm) of the static GNSS vector in the legend",
    ),
    map_limits: Tuple[float, float, float, float] = typer.Option(
        [None, None, None, None],
        "-ml",
        "--map-limits",
        help=(
            "Specify map limits [W,E,N,S] from edges of plotted features. "
            "eg: 0.5 0.5 0.5 0.5 gives a 0.5 degree buffer on each side. "
            "Negative numbers will cut off plotted features."
        ),
    ),
    max_slip: float = typer.Option(
        None,
        "-ms",
        "--max-slip",
        help="Set the maximum slip to be displayed on the plot",
    ),
    scale: float = typer.Option(
        None,
        "-sf",
        "--scale-factor",
        help="Scale factor for static GNSS vector lengths (larger means shorter vectors)",
    ),
):
    # set default data type
    chosen_data_types: List[str]
    chosen_data_types = [d.value for d in data_types]

    # Set and validate files
    paths_to_validate = []
    directory = directory.resolve()
    paths_to_validate += [directory]
    paths_to_validate += [gcmt_tensor_file]
    segments_file = directory / "segments_data.json"
    paths_to_validate += [segments_file]
    paths_to_validate += [directory / "Solucion.txt"]
    for dt in chosen_data_types:
        paths_to_validate += [directory / DEFAULT_MANAGEMENT_FILES[dt]]
    validate_files(paths_to_validate)

    # get default directories
    default_directories = default_dirs(config_path=config_file)

    # get the tensor information
    tensor_info = get_tensor(cmt_file=gcmt_tensor_file)

    # get the segments data and connections
    with open(segments_file) as sf:
        segments_data = json.load(sf)
    segments = segments_data["segments"]
    rise_time = segments_data["rise_time"]
    connections = None
    if "connections" in segments_data:
        connections = segments_data["connections"]
    point_sources = point_sources_param(
        segments, tensor_info, rise_time, connections=connections
    )

    # get the velocity model
    if not (directory / "velmodel_data.json").exists():
        vel_model = select_velmodel(
            tensor_info=tensor_info,
            default_dirs=default_directories,
            directory=directory,
        )
    else:
        with open(directory / "velmodel_data.json") as v:
            vel_model = json.load(v)

    # get the solution
    solution = read_solution_static_format(segments, data_dir=directory)

    # get files with trace information
    if ManagedDataTypes.gps in chosen_data_types:
        names, lats, lons, observed, synthetic, error = retrieve_gps(
            directory=directory
        )
        stations_gps = zip(names, lats, lons, observed, synthetic, error)
    else:
        stations_gps = None
    if ManagedDataTypes.cgps in chosen_data_types:
        with open(directory / DEFAULT_MANAGEMENT_FILES["cgps"]) as c:
            traces_info_cgps = json.load(c)
    else:
        traces_info_cgps = None
    if ManagedDataTypes.strong in chosen_data_types:
        with open(directory / DEFAULT_MANAGEMENT_FILES["strong"]) as s:
            traces_info = json.load(s)
    else:
        traces_info = None

    PlotMap_KML(
        tensor_info=tensor_info,
        segments=segments,
        point_sources=point_sources,
        solution=solution,
        default_dirs=default_directories,
        stations_str=traces_info,
        stations_gps=stations_gps,
        stations_cgps=traces_info_cgps,
        max_slip=max_slip,
        legend_len=legend_length,
        scale=scale,
        limits=list(map_limits),
        evID=event_id,
        directory=directory,
    )


@app.command(help="Create map plots")
def map(
    directory: pathlib.Path = typer.Argument(..., help="Path to the data directory"),
    gcmt_tensor_file: pathlib.Path = typer.Argument(
        ..., help="Path to the GCMT moment tensor file"
    ),
    config_file: pathlib.Path = typer.Option(
        CONFIG_PATH, "-c", "--config-file", help="Path to config file"
    ),
    data_types: List[ManagedDataTypes] = typer.Option(
        [],
        "-t",
        "--data-type",
        help="Data types to plot misfit for, default is []",
    ),
    legend_length: float = typer.Option(
        None,
        "-ll",
        "--legend-length",
        help="Set the length (cm) of the static GNSS vector in the legend",
    ),
    map_limits: Tuple[float, float, float, float] = typer.Option(
        [None, None, None, None],
        "-ml",
        "--map-limits",
        help=(
            "Specify map limits [W,E,N,S] from edges of plotted features. "
            "eg: 0.5 0.5 0.5 0.5 gives a 0.5 degree buffer on each side. "
            "Negative numbers will cut off plotted features."
        ),
    ),
    max_slip: float = typer.Option(
        None,
        "-ms",
        "--max-slip",
        help="Set the maximum slip to be displayed on the plot",
    ),
    scale: float = typer.Option(
        None,
        "-sf",
        "--scale-factor",
        help="Scale factor for static GNSS vector lengths (larger means shorter vectors)",
    ),
):
    # set default data type
    chosen_data_types: List[str]
    chosen_data_types = [d.value for d in data_types]

    # Set and validate files
    paths_to_validate = []
    directory = directory.resolve()
    paths_to_validate += [directory]
    paths_to_validate += [gcmt_tensor_file]
    segments_file = directory / "segments_data.json"
    paths_to_validate += [segments_file]
    paths_to_validate += [directory / "Solucion.txt"]
    for dt in chosen_data_types:
        paths_to_validate += [directory / DEFAULT_MANAGEMENT_FILES[dt]]
    validate_files(paths_to_validate)

    # get default directories
    default_directories = default_dirs(config_path=config_file)

    # get the tensor information
    tensor_info = get_tensor(cmt_file=gcmt_tensor_file)

    # get the segments data and connections
    with open(segments_file) as sf:
        segments_data = json.load(sf)
    segments = segments_data["segments"]
    rise_time = segments_data["rise_time"]
    connections = None
    if "connections" in segments_data:
        connections = segments_data["connections"]
    point_sources = point_sources_param(
        segments, tensor_info, rise_time, connections=connections
    )

    # get the solution
    solution = read_solution_static_format(segments, data_dir=directory)

    # get files with trace information
    if ManagedDataTypes.gps in chosen_data_types:
        names, lats, lons, observed, synthetic, error = retrieve_gps(
            directory=directory
        )
        stations_gps = zip(names, lats, lons, observed, synthetic, error)
    else:
        stations_gps = None
    if ManagedDataTypes.cgps in chosen_data_types:
        with open(directory / DEFAULT_MANAGEMENT_FILES["cgps"]) as c:
            traces_info_cgps = json.load(c)
    else:
        traces_info_cgps = None
    if ManagedDataTypes.strong in chosen_data_types:
        with open(directory / DEFAULT_MANAGEMENT_FILES["strong"]) as s:
            traces_info = json.load(s)
    else:
        traces_info = None

    PlotMap(
        tensor_info=tensor_info,
        segments=segments,
        point_sources=point_sources,
        solution=solution,
        default_dirs=default_directories,
        files_str=traces_info,
        stations_gps=stations_gps,
        stations_cgps=traces_info_cgps,
        max_slip=max_slip,
        legend_len=legend_length,
        scale=scale,
        limits=list(map_limits),
        directory=directory,
    )


@app.command(help="Create the plots used by the NEIC at the USGS")
def neic(
    directory: pathlib.Path = typer.Argument(..., help="Path to the data directory"),
    auto_size: bool = typer.Option(
        False,
        "-a",
        "--auto-size",
        help="Whether to automatically scale the plots",
    ),
    checkerboard: bool = typer.Option(
        False,
        "-cb",
        "--checkerboard",
        help="Plot comparisons for the checkerboard test",
    ),
    config_file: pathlib.Path = typer.Option(
        CONFIG_PATH, "-c", "--config-file", help="Path to config file"
    ),
    gcmt_tensor_file: pathlib.Path = typer.Option(
        None, "-g", "--gcmt", help="Path to the GCMT moment tensor file"
    ),
    data_types: List[ManagedDataTypes] = typer.Option(
        [],
        "-t",
        "--data-type",
        help="Data types to plot misfit for, default is []",
    ),
    downloads: bool = typer.Option(
        False,
        "-d",
        "--downloads",
        help="Create downloads to be displayed on the USGS event pages",
    ),
    event_id: str = typer.Option(
        None,
        "-e",
        "--event-id",
        help="The event's name/ID",
    ),
    ffm_solution: bool = typer.Option(
        False,
        "-ffms",
        "--ffm-solution",
        help="Plot the FFM solution slip maps and rise time",
    ),
    label: bool = typer.Option(
        False,
        "-l",
        "--label",
        help="Label local stations with their IDs on the map (for strong, cgps, gps data)",
    ),
    legend_length: float = typer.Option(
        None,
        "-ll",
        "--legend-length",
        help="Set the length (cm) of the static GNSS vector in the legend",
    ),
    map_limits: Tuple[float, float, float, float] = typer.Option(
        [None, None, None, None],
        "-ml",
        "--map-limits",
        help=(
            "Specify map limits [W,E,N,S] from edges of plotted features. "
            "eg: 0.5 0.5 0.5 0.5 gives a 0.5 degree buffer on each side. "
            "Negative numbers will cut off plotted features."
        ),
    ),
    max_slip: float = typer.Option(
        None,
        "-ms",
        "--max-slip",
        help="Set the maximum slip to be displayed on the plot",
    ),
    moment_rate_time: int = typer.Option(
        None,
        "-mrt",
        "--momentrate-time",
        help="Specify a cutoff time for the Moment Rate plot",
    ),
    polygon: bool = typer.Option(
        False,
        "-p",
        "--polygon",
        help="Create a polygon to be used by ShakeMap",
    ),
    publish: bool = typer.Option(
        False,
        "-pub",
        "--publish",
        help="Make files for publishing the finite fault product to the USGS event pages",
    ),
    scale: float = typer.Option(
        None,
        "-sf",
        "--scale-factor",
        help="Scale factor for static GNSS vector lengths (larger means shorter vectors)",
    ),
    separate_planes: bool = typer.Option(
        False,
        "-sp",
        "--separate-planes",
        help="Include STFs for each plane separately in the moment rate plot",
    ),
    tensor: bool = typer.Option(
        False,
        "--tensor",
        help="Calculate the cumulative moment tensor of the finite fault model",
    ),
):
    # set default data type
    chosen_data_types: List[str]
    chosen_data_types = [d.value for d in data_types]

    # Set and validate files
    paths_to_validate = []
    directory = directory.resolve()
    paths_to_validate += [directory]
    tensor_file = directory / "tensor_info.json"
    paths_to_validate += [tensor_file]
    segments_file = directory / "segments_data.json"
    paths_to_validate += [segments_file]
    paths_to_validate += [directory / "Solucion.txt"]
    for dt in chosen_data_types:
        paths_to_validate += [directory / DEFAULT_MANAGEMENT_FILES[dt]]
    validate_files(paths_to_validate)

    # get default directories
    default_directories = default_dirs(config_path=config_file)

    # get the tensor information
    with open(tensor_file) as tf:
        tensor_info = json.load(tf)

    # get the segments data and connections
    with open(segments_file) as sf:
        segments_data = json.load(sf)
    segments = segments_data["segments"]
    rise_time = segments_data["rise_time"]
    connections = None
    if "connections" in segments_data:
        connections = segments_data["connections"]
    point_sources = point_sources_param(
        segments, tensor_info, rise_time, connections=connections
    )

    # get the velocity model
    if not (directory / "velmodel_data.json").exists():
        vel_model = select_velmodel(
            tensor_info=tensor_info,
            default_dirs=default_directories,
            directory=directory,
        )
    else:
        with open(directory / "velmodel_data.json") as v:
            vel_model = json.load(v)

    # get the solution
    solution = read_solution_static_format(segments, data_dir=directory)

    # get files with trace information
    if ManagedDataTypes.gps in chosen_data_types:
        names, lats, lons, observed, synthetic, error = retrieve_gps(
            directory=directory
        )
        stations_gps = zip(names, lats, lons, observed, synthetic, error)
    else:
        stations_gps = None
    if ManagedDataTypes.cgps in chosen_data_types:
        with open(directory / DEFAULT_MANAGEMENT_FILES["cgps"]) as c:
            traces_info_cgps = json.load(c)
    else:
        traces_info_cgps = None
    if ManagedDataTypes.strong in chosen_data_types:
        with open(directory / DEFAULT_MANAGEMENT_FILES["strong"]) as s:
            traces_info = json.load(s)
    else:
        traces_info = None

    # generate other formats of the solution
    static_to_fsp(
        tensor_info=tensor_info,
        segments_data=segments_data,
        used_data=chosen_data_types,
        vel_model=vel_model,
        solution=solution,
        directory=directory,
    )
    static_to_srf(
        tensor_info=tensor_info,
        segments_data=segments_data,
        used_data=chosen_data_types,
        vel_model=vel_model,
        solution=solution,
        directory=directory,
    )

    # performing specific requests
    if checkerboard:
        input_model = load_ffm_model(
            segments_data=segments_data,
            point_sources=point_sources,
            option="fault&rise_time.txt",
            directory=directory,
        )
        PlotComparisonMap(
            tensor_info=tensor_info,
            segments=segments_data,
            point_sources=point_sources,
            input_model=input_model,
            solution=solution,
            max_val=max_slip,
            directory=directory,
        )
        PlotSlipDist_Compare(
            segments=segments_data,
            point_sources=point_sources,
            input_model=input_model,
            solution=solution,
            max_val=max_slip,
            directory=directory,
        )
    if tensor:
        calculate_cumulative_moment_tensor(solution=solution, directory=directory)
    if downloads:
        if gcmt_tensor_file:
            write_CMTSOLUTION_file(pdefile=gcmt_tensor_file, directory=directory)
        else:
            logging.warn(
                "No gcmt_tensor_file specified. Skipping writing CMTSOLUTION file."
            )
        write_Coulomb_file(directory=directory, eventID=event_id)
        write_Okada_displacements(directory=directory)
        make_waveproperties_json(directory=directory)
    if ffm_solution:
        shear = shear_modulous(point_sources, velmodel=vel_model)
        plot_ffm_sol(
            tensor_info=tensor_info,
            segments_data=segments_data,
            point_sources=point_sources,
            shear=shear,
            solution=solution,
            default_dirs=default_directories,
            autosize=auto_size,
            mr_time=moment_rate_time,
            files_str=traces_info,
            stations_gps=stations_gps,
            stations_cgps=traces_info_cgps,
            max_val=max_slip,
            legend_len=legend_length,
            scale=scale,
            limits=list(map_limits),
            separate_planes=separate_planes,
            label_stations=label,
            directory=directory,
        )
    if ManagedDataTypes.insar in chosen_data_types:
        insar_data = get_insar(data_dir=directory)
        if "ascending" in insar_data:
            for scene in range(len(insar_data["ascending"])):
                insar_points = insar_data["ascending"][scene]["points"]
                PlotInsar(
                    tensor_info=tensor_info,
                    segments=segments,
                    point_sources=point_sources,
                    solution=solution,
                    insar_points=insar_points,
                    scene=str(scene),
                    los="ascending",
                    directory=directory,
                )
        if "descending" in insar_data:
            for scene in range(len(insar_data["descending"])):
                insar_points = insar_data["descending"][scene]["points"]
                PlotInsar(
                    tensor_info=tensor_info,
                    segments=segments,
                    point_sources=point_sources,
                    solution=solution,
                    insar_points=insar_points,
                    scene=str(scene),
                    los="descending",
                    directory=directory,
                )
    if polygon:
        shakemap_polygon(
            segments=segments,
            point_sources=point_sources,
            solution=solution,
            tensor_info=tensor_info,
            evID=event_id,
            directory=directory,
        )

    # plot misfit
    plot_misfit(used_data_type=chosen_data_types, directory=directory)

    # set up plot directory
    if not (directory / "plots").exists():
        os.mkdir(directory / "plots")
    plot_files = glob(str(directory) + "/*png")
    for plot_file in plot_files:
        basename = os.path.basename(plot_file)
        move(
            pathlib.Path(plot_file).resolve(),
            (directory / "plots" / basename).resolve(),
        )

    # set up publishing directory
    if publish:
        if event_id is None:
            raise ValueError("event_id must be specified in order to run publish")
        temporary_file_reorganization_for_publishing(evID=event_id, directory=directory)
