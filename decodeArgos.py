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
    df = df.set_index("time")
    ds = df.to_xarray()
    ds.to_netcdf(fnNC)

def procFile(fn:str) -> pd.DataFrame:
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
            matches = re.fullmatch(
                    r"(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(.)\s+(.)\s+(\d{4})-(\d{2})-(\d{2})\s+(\d{2}):(\d{2}):(\d{2})\s+(\d+[.]\d+)\s+(\d+[.]\d+)\s+(\d+[.]\d+)\s+(\d+)", 
                               line)
            if not matches: continue
            records["ident"].append(int(matches[2]))
            records["nLines"].append(int(matches[3]))
            records["nBytes"].append(int(matches[4]))
            records["satellite"].append(matches[5])
            records["locationClass"].append(matches[6])
            records["time"].append(datetime.datetime(
                int(matches[7]), 
                int(matches[8]), 
                int(matches[9]),
                int(matches[10]),
                int(matches[11]),
                int(matches[12])))
            records["lat"].append(float(matches[13]))
            records["lon"].append(float(matches[14]))
            records["altitude"].append(float(matches[15]))
            records["frequency"].append(float(matches[16]))

    if len(records["ident"]) == 0: return None
    df = pd.DataFrame(records)
    return df

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("filename", type=str, nargs="+", help="ARGOS messages to decode")
    parser.add_argument("--nc", type=str, default="tpw.nc", help="Output NetCDF filename")
    args = parser.parse_args()

    processFiles(args.nc, args.filename)
