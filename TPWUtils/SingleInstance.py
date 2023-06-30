#! /usr/bin/env python3
#
# Using socket listener to check there is only one listener for a specified port
#
# June-2022, Pat Welch, pat@mousebrains.com

import socket
import logging
import os
import sys

class SingleInstance: # Must be used with with statement
    def __init__(self, key:str = None) -> None:
        self.__key = os.path.abspath(os.path.expanduser(sys.argv[0])) if key is None else key

    def __enter__(self):
        try:
            s = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
            s.bind('\0' + self.__key) # Abstract socket by prefixing with a null
            self.__socket = s
        except:
            logging.exception("Unable to connect to %s", self.__key)
            self.__socket = None
            sys.exit(0)

    def __exit__(self, excType, excValue, excTraceback) -> None:
        self.__socket = None

if __name__ == "__main__":
    from argparse import ArgumentParser
    import time

    parser = ArgumentParser()
    parser.add_argument("--uniqueName", type=str,
            help="Single instance unique keyword for locking a process")
    parser.add_argument("--dt", type=float, default=100, help="Time to sleep")
    args = parser.parse_args()
    
    logging.basicConfig(level=logging.DEBUG, format="%(asctime)s %(levelname)s: %(message)s")
   
    try:
        with SingleInstance(args.uniqueName) as single:
            logging.info("Sleeping for %s seconds", args.dt)
            time.sleep(args.dt)
            logging.info("Done sleeping")
    except SystemExit:
        logging.error("Single instance failure")
    except:
        logging.exception("Unexpected exception")
