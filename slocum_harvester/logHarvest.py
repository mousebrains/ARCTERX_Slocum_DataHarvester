#! /usr/bin/env python3
#
# Harvest records from the log files and save to a NetCDF for further analysis
#
# June-2023, Pat Welch, pat@mousebrains.com

import glob
import logging
import xarray as xr
import pandas as pd
import numpy as np
import os.path
import datetime
import re
import sys
from slocum_harvester.slocum_utils import mkDegrees_scalar as mkDegrees

def parseLogFile(fn, glider) -> pd.DataFrame:
    logging.debug("Loading %s %s", fn, glider)
    info = {}
    times = []
    with open(fn, "rb") as fp:
        currTime = None
        missionTime = None
        for line in fp:
            line = line.strip()
            try:
                line = str(line, "utf-8")
            except UnicodeDecodeError:
                continue
            matches = re.match(r"^Vehicle Name:\s+(\w+)$", line)
            if matches: glider = matches[1]

            matches = re.match(r"^Curr Time:\s+\w+\s+(\w+\s+\d+\s+\d+:\d+:\d+\s+\d+)\s+MT:\s+(\d+)", line)
            if matches:
                try:
                    currTime = datetime.datetime.strptime(matches[1], "%b %d %H:%M:%S %Y") \
                            .replace(tzinfo=datetime.timezone.utc).timestamp()
                except ValueError:
                    logging.warning("Invalid timestamp in %s: %s", fn, matches[1])
                    continue
                missionTime = int(matches[2])
                continue
            matches = re.match(r"^GPS Location:\s+(\d+[.]\d+)\s+N\s+(\d+[.]\d+)\s+E\s+measured\s+([+-]?\d+[.]?\d*[e]?[+-]?\d*) secs ago$", line)
            if matches:
                if currTime is None: continue
                lat = mkDegrees(float(matches[1]))
                lon = mkDegrees(float(matches[2]))
                dt = float(matches[3])
                if dt > 1e300 or abs(lat) > 90 or abs(lon) > 180: continue
                key="GPS"
                if key not in info: info[key] = []
                t = currTime - dt
                times.append(t)
                info[key].append((t, lat, lon))
                continue
            matches = re.match(r"^sensor:(\w+)[(](.+)[)]=([+-]?\d+[.]?\d*[e]?[+-]?\d*)\s+(\d+[.]?\d*[e]?[+-]?\d*) secs ago$", line)
            if matches:
                if currTime is None: continue
                name = matches[1]
                val = float(matches[3])
                if name.endswith("_lat") or name.endswith("_lon"): val = mkDegrees(val)
                dt = float(matches[4])
                if dt > 1e300: continue
                if name not in info: info[name] = []
                t = currTime - dt
                times.append(t)
                info[name].append((t, val))

    df = pd.DataFrame()
    df["t"] = np.unique(np.round(np.array(times), -2))
    df["glider"] = np.array([glider] * df.t.size)
    for key in sorted(info):
        if key == "GPS":
            df["lat"] = np.empty(df.t.shape) * np.nan
            df["lon"] = np.empty(df.t.shape) * np.nan
        else:
            df[key] = np.empty(df.t.shape) * np.nan

    for key in sorted(info):
        for row in info[key]:
            row = np.array(row)
            index = np.argmin(abs(df.t - row[0].round(-2)))
            if key == "GPS":
                df.loc[index, "lat"] = row[1]
                df.loc[index, "lon"] = row[2]
            else:
                df.loc[index, key] = row[1]
    df["t"] = df["t"].astype("datetime64[s]")
    return df

def processFiles(filenames:list, t0:str, nc:str):
    items = []
    for fn in sorted(filenames):
        fields = os.path.basename(fn).split("_")
        if len(fields) < 2: continue
        if t0 and fields[1] < t0: continue
        glider = fields[0]
        a = parseLogFile(fn, glider)
        if a is not None and a.t.size:
            items.append(a)

    if not items:
        logging.warning("No valid log files found, writing empty %s", nc)
        xr.Dataset().to_netcdf(nc)
        return

    ds = xr.Dataset.from_dataframe(pd.concat(items, ignore_index=True))
    logging.info("Writing %s to %s", ds.sizes, nc)
    ds.to_netcdf(nc)

def main():
    from argparse import ArgumentParser
    from TPWUtils import Logger

    parser = ArgumentParser(description="Harvest Slocum glider log files into NetCDF")
    Logger.addArgs(parser)
    parser.add_argument("--t0", type=str, default="20230527T000000", help="Earliest time to consider")
    parser.add_argument("--nc", type=str, default="log.nc", help="Output NetCDF filename")
    parser.add_argument("filename", type=str, nargs="*", help="Log files to parse")
    args = parser.parse_args()

    Logger.mkLogger(args, fmt="%(asctime)s %(levelname)s: %(message)s")

    if not args.filename:
        args.filename = glob.glob("osu68?/logs/*.log")

    processFiles(args.filename, args.t0, args.nc)

if __name__ == "__main__":
    main()
