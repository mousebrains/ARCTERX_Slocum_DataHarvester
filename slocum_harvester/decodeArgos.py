#! /usr/bin/env python3
#
# Work on decoding the hex ARGOS message

import logging
import re
import pandas as pd
import datetime

# ARGOS message line format:
# Fields: unknown, ident, nLines, nBytes, satellite, locationClass,
#         YYYY-MM-DD HH:MM:SS, lat, lon, altitude, frequency
_ARGOS_RE = re.compile(r"""
    (\d+)       \s+   # field 1 (unused)
    (\d+)       \s+   # ident
    (\d+)       \s+   # nLines
    (\d+)       \s+   # nBytes
    (.)         \s+   # satellite ID
    (.)         \s+   # location class
    (\d{4}) - (\d{2}) - (\d{2}) \s+  # date YYYY-MM-DD
    (\d{2}) : (\d{2}) : (\d{2}) \s+  # time HH:MM:SS
    (\d+[.]\d+) \s+   # latitude
    (\d+[.]\d+) \s+   # longitude
    (\d+[.]\d+) \s+   # altitude
    (\d+)              # frequency
""", re.VERBOSE)

def processFiles(fnNC:str, filenames:list):
    records = []
    for fn in filenames:
        df = procFile(fn)
        if df is not None: records.append(df)

    if not records:
        logging.warning("No ARGOS records found in %s files, writing empty %s",
                        len(filenames), fnNC)
        pd.DataFrame().to_xarray().to_netcdf(fnNC)
        return

    df = pd.concat(records, ignore_index=True)
    df = df.set_index("time")
    ds = df.to_xarray()
    ds.to_netcdf(fnNC)

def procFile(fn:str) -> pd.DataFrame | None:
    records = dict(
            ident = [],
            nLines = [],
            nBytes = [],
            satellite = [],
            locationClass = [],
            time = [],
            lat = [],
            lon = [],
            altitude = [],
            frequency = [],
            )

    with open(fn, "r") as fp:
        for line in fp:
            line = line.strip()
            matches = _ARGOS_RE.fullmatch(line)
            if not matches: continue
            try:
                t = datetime.datetime(
                    int(matches[7]), int(matches[8]), int(matches[9]),
                    int(matches[10]), int(matches[11]), int(matches[12]),
                    tzinfo=datetime.timezone.utc)
            except ValueError:
                continue
            records["ident"].append(int(matches[2]))
            records["nLines"].append(int(matches[3]))
            records["nBytes"].append(int(matches[4]))
            records["satellite"].append(matches[5])
            records["locationClass"].append(matches[6])
            records["time"].append(t)
            records["lat"].append(float(matches[13]))
            records["lon"].append(float(matches[14]))
            records["altitude"].append(float(matches[15]))
            records["frequency"].append(float(matches[16]))

    if len(records["ident"]) == 0: return None
    df = pd.DataFrame(records)
    return df

def main():
    from argparse import ArgumentParser
    parser = ArgumentParser(description="Decode ARGOS satellite messages into NetCDF")
    parser.add_argument("filename", type=str, nargs="+", help="ARGOS messages to decode")
    parser.add_argument("--nc", type=str, default="tpw.nc", help="Output NetCDF filename")
    args = parser.parse_args()

    processFiles(args.nc, args.filename)

if __name__ == "__main__":
    main()
