#! /usr/bin/env python3
#
# Update a fresh copy of the glider files
# Harvest the log files into NetCDF files
# Harvest the s?d and t?d files into NetCDF files
# Combine the harvested information into a NetCDF file for distribution
#
# June-2023, Pat Welch, pat@mousebrains.com

from argparse import ArgumentParser
import subprocess
from TPWUtils import Logger
import logHarvest
import decodeArgos
from mkCombined import mkCombo
import logging
import glob
import os
import shutil
import sys

def myMkDir(directory:str) -> None:
    if os.path.isdir(directory): return
    logging.info("Creating %s", directory)
    os.makedirs(directory, mode=0o755, exist_ok=True) # exist_ok for race conditions

def runCmd(cmd, returnCodes:tuple=(0,)) -> bool:
    logging.debug("cmd %s", cmd)
    a = subprocess.run(cmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    if a.returncode in returnCodes:
        if a.stdout:
            logging.debug("Ran %s", " ".join(cmd))
            try:
                logging.debug("%s", str(a.stdout, "UTF-8"))
            except Exception:
                logging.debug("%s", a.stdout)
        return True

    if a.stdout:
        logging.warning("Failed (rc=%s) %s", a.returncode, " ".join(cmd))
        try:
            logging.warning("%s", str(a.stdout, "UTF-8"))
        except Exception:
            logging.warning("%s", a.stdout)
    return False

def main():
    parser = ArgumentParser()
    Logger.addArgs(parser)
    parser.add_argument("--hostname", type=str, default="fs2", help="Hostname to rsync from")
    parser.add_argument("--src", type=str, action="append", help="Input source for rsync")
    parser.add_argument("--rsync", type=str, default="/usr/bin/rsync", help="rsync path")
    parser.add_argument("--dbd2netCDF", type=str, default="/usr/local/bin/dbd2netCDF",
                        help="dbd2netCDF path")
    parser.add_argument("--glider", type=str, action="append", help="Glider(s) to process")
    parser.add_argument("--t0", type=str, default="20230527T000000",
                        help="Earliest time to consider in logfiles")
    parser.add_argument("--tgt", type=str, default=".", help="Output directory for rsync")
    parser.add_argument("--output", type=str, default=".",
                        help="Output directory for log/sbd/tbd NetCDF files")
    parser.add_argument("--combo", type=str, default=".",
                        help="Output directory for combined NetCDF files")
    args = parser.parse_args()

    if args.glider is None:
        args.glider = ("osu684", "osu685", "osu686")

    if args.src is None:
        args.src = (
                f"{args.hostname}:/data/Dockserver/gliderfmc0/osu68?",
                f"{args.hostname}:/data/tpw/cache",
                )

    args.rsync = shutil.which(args.rsync) or shutil.which("rsync")
    args.dbd2netCDF = shutil.which(args.dbd2netCDF) or shutil.which("dbd2netCDF")
    args.tgt = os.path.abspath(os.path.expanduser(args.tgt))
    args.output = os.path.abspath(os.path.expanduser(args.output))
    args.combo = os.path.abspath(os.path.expanduser(args.combo))

    Logger.mkLogger(args, fmt="%(asctime)s %(levelname)s: %(message)s")

    if not args.rsync:
        logging.error("rsync not found")
        sys.exit(2)
    if not args.dbd2netCDF:
        logging.error("dbd2netCDF not found")
        sys.exit(2)

    try:
        myMkDir(args.tgt)
        myMkDir(args.output)
        myMkDir(args.combo)

        cmdBase = [args.rsync, "--compress", "--archive", "--delete"]
        if args.verbose or args.debug: cmdBase.append("--verbose")
        for src in args.src:
            cmd = list(cmdBase) # Make a local copy
            cmd.extend([src, args.tgt])
            if not runCmd(cmd): sys.exit(1)

        # Grab ARGOS fixes
        cmd = list(cmdBase)
        cmd.extend(["--exclude", "ornitela", "ceoas:/home/server/ftp/pub/incoming", args.tgt])
        if not runCmd(cmd): sys.exit(1)

        decodeArgos.processFiles(
                os.path.join(args.output, "argos.nc"),
                glob.glob(os.path.join(args.tgt, "incoming/*.TXT")))

        logHarvest.processFiles(
                glob.glob(os.path.join(args.tgt, "*/logs/*.log")),
                args.t0,
                os.path.join(args.output, "log.nc"))

        cmdBase = [args.dbd2netCDF, "--repair", "--skipFirst",
                   "--cache", os.path.join(args.tgt, "cache")]
        if args.verbose or args.debug: cmdBase.append("--verbose")
        cmdBase.append("--output")

        for glider in args.glider:
            log = os.path.join(args.output, "log.nc")
            flt = os.path.join(args.output, f"flt.{glider}.nc")
            sci = os.path.join(args.output, f"sci.{glider}.nc")
            logging.info("Working on %s", glider)
            cmd = list(cmdBase)
            cmd.append(flt)
            cmd.extend(glob.glob(os.path.join(args.tgt, glider, "from-glider", "*.s?d")))
            if not runCmd(cmd): sys.exit(3)
            cmd = list(cmdBase)
            cmd.append(sci)
            cmd.extend(glob.glob(os.path.join(args.tgt, glider, "from-glider", "*.t?d")))
            if not runCmd(cmd): sys.exit(4)

            if not mkCombo(glider,
                    os.path.join(args.combo, f"{glider}.nc"),
                    log, flt, sci):
                sys.exit(5)
    except Exception:
        logging.exception("Unexpected exception, %s", args)
        sys.exit(1)

if __name__ == "__main__":
    main()
