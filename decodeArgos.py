#! /usr/bin/env python3
#
# Work on decoding the hex ARGOS message

import re
import numpy as np
import pandas as pd
import datetime

def processFiles(fnNC:str, filenames:list):
    records = []
    for fn in filenames:
        df = procFile(fn)
        if df is not None: records.append(df)

    df = pd.concat(records, ignore_index=True)
    df = df.set_index("t")
    ds = df.to_xarray()
    ds.to_netcdf(fnNC)

def procFile(fn:str) -> pd.DataFrame:
    ident = []
    time = []
    lat = []
    lon = []

    with open(fn, "r") as fp:
        for line in fp:
            line = line.strip()
            matches = re.fullmatch(
                    r"(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(.)\s+(.)\s+(\d{4})-(\d{2})-(\d{2})\s+(\d{2}):(\d{2}):(\d{2})\s+(\d+[.]\d+)\s+(\d+[.]\d+)\s+(\d+[.]\d+)\s+(\d+)", 
                               line)
            if not matches: continue
            ident.append(int(matches[2]))
            time.append(datetime.datetime(
                int(matches[7]), 
                int(matches[8]), 
                int(matches[9]),
                int(matches[10]),
                int(matches[11]),
                int(matches[12])))
            lat.append(float(matches[13]))
            lon.append(float(matches[14]))

    if len(ident) == 0: return None
    df = pd.DataFrame()
    df["ident"] = ident
    df["t"] = time
    df["lat"] = lat
    df["lon"] = lon
    return df

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("filename", type=str, nargs="+", help="ARGOS messages to decode")
    parser.add_argument("--nc", type=str, default="tpw.nc", help="Output NetCDF filename")
    args = parser.parse_args()

    processFiles(args.nc, args.filename)
