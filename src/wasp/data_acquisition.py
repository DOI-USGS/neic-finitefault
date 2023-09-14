# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 17:34:16 2017

@author: pk
"""

# Manual data acquisition


import os
import pathlib
import queue
import threading
import time
from typing import Dict, List, Literal, Union

from obspy.clients.fdsn import Client  # type: ignore
from obspy.clients.iris import Client as IRIS_Client  # type: ignore
from obspy.core.inventory.channel import Channel  # type: ignore
from obspy.core.inventory.inventory import Inventory  # type: ignore
from obspy.core.inventory.util import Latitude, Longitude  # type: ignore
from obspy.core.utcdatetime import UTCDateTime  # type: ignore

import wasp.seismic_tensor as tensor


def acquisition(
    event_time: UTCDateTime,
    lat_ep: float,
    lon_ep: float,
    depth: float,
    data_to_use: List[Literal["strong", "tele"]],
    minutes_before: int = 3,
    minutes_after: int = 67,
    strong_before: int = 1,
    strong_after: int = 5,
    tele_before: int = 10,
    tele_after: int = 120,
    waveform_directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Get station data from IRIS

    :param event_time: The time of the event
    :type event_time: UTCDateTime
    :param lat_ep: The latitude of the epicenter
    :type lat_ep: float
    :param lon_ep: The longitude of the epicenter
    :type lon_ep: float
    :param depth: The depth of the epicenter
    :type depth: float
    :param data_to_use: The type of data to use
    :type data_to_use: List[str]
    :param minutes_before: Minutes pre-event included in waveform query, defaults to 3
    :type minutes_before: int, optional
    :param minutes_after: Minutes post-event included in waveform query, defaults to 67
    :type minutes_after: int, optional
    :param strong_before: Minutes pre-event included in station query, defaults to 1
    :type strong_before: int, optional
    :param strong_after: Minutes post-event included in station query, defaults to 5
    :type strong_after: int, optional
    :param tele_before: Minutes pre-event included in station query, defaults to 10
    :type tele_before: int, optional
    :param tele_after: Minutes post-event included in station query, defaults to 120
    :type tele_after: int, optional
    :param waveform_directory: _description_, defaults to "."
    :type waveform_directory: Union[pathlib.Path, str], optional
    """
    t1 = event_time - minutes_before * 60
    t2 = event_time + minutes_after * 60
    client_iris = Client("IRIS")
    try:
        client_gfz = Client("GFZ")  # For stations belonging to CX network
    except:
        pass
    inventory = Inventory([], None)
    if "strong" in data_to_use:
        networks = "C,C1,II,IU"
        try:
            inventory = client_iris.get_stations(
                starttime=event_time - strong_before * 60,
                endtime=event_time + strong_after * 60,
                network=networks,
                channel="HN*",
                level="response",
                maxradius=10,
                latitude=lat_ep,
                longitude=lon_ep,
            )
        except Exception as e:
            print(e)
        networks = "CX"
        try:
            inventory = inventory + client_gfz.get_stations(
                starttime=event_time - strong_before * 60,
                endtime=event_time + strong_after * 60,
                network=networks,
                channel="HL*",
                level="response",
                maxradius=10,
                latitude=lat_ep,
                longitude=lon_ep,
            )
            print("we have the inventory")
        except:
            pass
    inventory_tele = Inventory([], None)
    if "tele" in data_to_use:
        networks = "II,G,IU,GE"
        max_dist = 90
        inventory_tele = client_iris.get_stations(
            starttime=event_time - tele_before * 60,
            endtime=event_time + tele_after * 60,
            network=networks,
            channel="BH*",
            level="response",
            minradius=30,
            maxradius=max_dist,
            latitude=lat_ep,
            longitude=lon_ep,
        )
    cola: queue.Queue = queue.Queue()
    for i in range(10):
        worker = threading.Thread(target=wrapper, args=(cola,), daemon=True)
        worker.start()

    iris_client = IRIS_Client()
    for network in inventory_tele:
        netwk = network.code
        for station in network:
            statn = station.code
            for canal in station:
                loc_code = canal.location_code
                channel = canal.code
                sac_dict = __get_channel_information_manual(
                    canal, lat_ep, lon_ep, depth
                )
                cola.put(
                    [
                        client_iris,
                        iris_client,
                        netwk,
                        statn,
                        loc_code,
                        canal,
                        sac_dict,
                        t1,
                        t2,
                        "teleseismic",
                        waveform_directory,
                    ]
                )

    for network in inventory:
        netwk = network.code
        for station in network:
            statn = station.code
            if statn in ["AC02"]:
                continue
            for canal in station:
                loc_code = canal.location_code
                channel = canal.code
                sac_dict = __get_channel_information_manual(
                    canal, lat_ep, lon_ep, depth
                )
                new_client = client_iris if not netwk == "CX" else client_gfz
                cola.put(
                    [
                        new_client,
                        iris_client,
                        netwk,
                        statn,
                        loc_code,
                        canal,
                        sac_dict,
                        t1,
                        t2,
                        "strong_motion",
                        waveform_directory,
                    ]
                )
    cola.join()
    time.sleep(5)
    return


def wrapper(q: queue.Queue):
    """Wrapper used to determine which worker should be used

    :param q: All of the arguments passed to the queue
    :type q: queue.Queue
    """
    while True:
        try:
            (
                client1,
                client2,
                netwk,
                statn,
                loc_code,
                canal,
                sac_dict,
                t1,
                t2,
                data_type,
                waveform_directory,
            ) = q.get(
                timeout=3
            )  # or whatever
        except Exception as e:
            print(f"Exception getting queue element: {e}")
            return
        channel = canal.code
        worker2(
            client1,
            netwk,
            statn,
            loc_code,
            channel,
            sac_dict,
            t1,
            t2,
            waveform_directory,
        )
        if data_type == "teleseismic":
            worker3(
                client2, netwk, statn, loc_code, channel, t1, t2, waveform_directory
            )
        else:
            worker4(canal, netwk, statn, loc_code, channel, t1, t2, waveform_directory)
        q.task_done()


def worker2(
    client: Union[Client, IRIS_Client],
    netwk: str,
    statn: str,
    loc_code: str,
    channel: str,
    sac_dict: Dict[str, Union[str, float]],
    time0: UTCDateTime,
    time1: UTCDateTime,
    waveform_directory: pathlib.Path,
):
    """Worker to get data from waveforms in sac format

    :param client: The obspy client
    :type client: Union[Client, IRIS_Client]
    :param netwk: The network code
    :type netwk: str
    :param statn: The station code
    :type statn: str
    :param loc_code: The location code
    :type loc_code: str
    :param channel: The channel code
    :type channel: str
    :param sac_dict: The sac dictionary
    :type sac_dict: Dict[str, Union[str, float]]
    :param time0: The start time of the query
    :type time0: UTCDateTime
    :param time1: The end time of the query
    :type time1: UTCDateTime
    :param waveform_directory: The directory where waveforms should be saved
    :type waveform_directory: pathlib.Path
    """
    try:
        st = client.get_waveforms(netwk, statn, loc_code, channel, time0, time1)
        name = "{}{}_{}{}.sac".format(netwk, statn, channel, loc_code)
        if channel[:2] in ["HL", "HN"] and len(st) > 1:
            return
        st[0].stats.sac = sac_dict
        st.write(str(waveform_directory / name), format="SAC", byteorder=0)
    except Exception as e:
        print(
            f"Exception in worker 2 for {netwk}, {statn}, {loc_code}, {channel}: "
            f"{e}"
        )
    return


def worker3(
    iris_client: IRIS_Client,
    netwk: str,
    statn: str,
    loc_code: str,
    channel: str,
    t1: UTCDateTime,
    t2: UTCDateTime,
    waveform_directory: pathlib.Path,
):
    """Worker to get instrument response information for teleseismic data

    :param client: The obspy IRIS client
    :type client: IRIS_Client
    :param netwk: The network code
    :type netwk: str
    :param statn: The station code
    :type statn: str
    :param loc_code: The location code
    :type loc_code: str
    :param channel: The channel code
    :type channel: str
    :param t1: The start time of the query
    :type t1: UTCDateTime
    :param t2: The end time of the query
    :type t2: UTCDateTime
    :param waveform_directory: The directory where waveforms should be saved
    :type waveform_directory: pathlib.Path
    """
    try:
        response_data = "SAC_PZs_{}_{}_{}_{}".format(netwk, statn, channel, loc_code)
        if len(loc_code) is 0:
            response_data = "SAC_PZs_{}_{}_{}___".format(netwk, statn, channel)
        iris_client.sacpz(
            netwk,
            statn,
            loc_code,
            channel,
            t1,
            t2,
            filename=str(waveform_directory / response_data),
        )
    except Exception as e:
        print(
            f"Exception in worker 3 for {netwk}, {statn}, {loc_code}, {channel}: "
            f"{e}"
        )
    return


def worker4(
    canal: Channel,
    netwk: str,
    statn: str,
    loc_code: str,
    channel: str,
    t1: UTCDateTime,
    t2: UTCDateTime,
    waveform_directory: pathlib.Path,
):
    """Worker to get instrument response information for strong motion data

    :param canal: The channel
    :type canal: Channel
    :param netwk: The network code
    :type netwk: str
    :param statn: The station code
    :type statn: str
    :param loc_code: The location code
    :type loc_code: str
    :param channel: The channel code
    :type channel: str
    :param t1: The start time of the query
    :type t1: UTCDateTime
    :param t2: The end time of the query
    :type t2: UTCDateTime
    :param waveform_directory: The directory where waveforms should be saved
    :type waveform_directory: pathlib.Path
    """
    try:
        response = canal.response
        response_data = "SAC_PZs_{}_{}_{}_{}".format(netwk, statn, channel, loc_code)
        if len(loc_code) is 0:
            response_data = "SAC_PZs_{}_{}_{}___".format(netwk, statn, channel)
        sacpz = response.get_sacpz()
        with open(waveform_directory / response_data, "w") as file:
            file.write("{}\n".format(sacpz))
    except Exception as e:
        print(
            f"Exception in worker 4 for {netwk}, {statn}, {loc_code}, {channel}: "
            f"{e}"
        )
    return


def __get_channel_information_manual(
    channel: Channel, lat_ep: float, lon_ep: float, depth: float
) -> Dict[str, Union[Latitude, Longitude, float, str]]:
    """Creates a dictionary with information and one of the
       components/channels of the station

    :param channel: The component/channel
    :type channel: Channel
    :param lat_ep: The epicenter latitude
    :type lat_ep: float
    :param lon_ep: The epicenter longitude
    :type lon_ep: float
    :param depth: The epicenter depth
    :type depth: float
    :return: The component/channel information dictionary
    :rtype: Dict[str, Union[Latitude, Longitude, float, str]]
    """
    station_lat = channel.latitude
    station_lon = channel.longitude
    loc_code = channel.location_code
    cmpaz = channel.azimuth
    cmpinc = 90 + channel.dip
    sac_dict = {
        "stla": station_lat,
        "stlo": station_lon,
        "khole": loc_code,
        "cmpaz": cmpaz,
        "cmpinc": cmpinc,
        "evla": lat_ep,
        "evlo": lon_ep,
        "evdp": depth,
    }
    return sac_dict
