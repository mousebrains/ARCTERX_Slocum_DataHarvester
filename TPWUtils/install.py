#! /usr/bin/env python3
#
# Install the service(s) and timer(s)
#
# Feb-2022, Pat Welch, pat@mousebrains.com
# June-2023, Pat Welch, pat@mousebrains.com updated for both root and user

from argparse import ArgumentParser
import logging
import subprocess
import os
import sys

def makeDirectory(dirname:str, args:ArgumentParser, qUser:bool=False) -> str:
    dirname = os.path.abspath(os.path.expanduser(dirname))
    if os.path.isdir(dirname): return dirname
    cmd = []
    if not qUser and not args.user: cmd.append(args.sudo)
    cmd.extend((args.mkdir, "-p", dirname))
    logging.info("Creating %s", " ".join(cmd))
    if not args.dryrun: subprocess.run(cmd, shell=False, check=True)
    return dirname

def stripComments(fn:str) -> str:
    lines = []
    with open(fn, "r") as fp:
        for line in fp.readlines():
            index = line.find("#")
            if index >= 0:
                line = line[:index]
            line = line.strip()
            if line: lines.append(line)
    return "\n".join(lines)

def needsToBeCopied(src:str, args:ArgumentParser) -> str:
    src = os.path.abspath(os.path.expanduser(src))
    tgt = os.path.join(args.serviceDirectory, os.path.basename(src))

    if not args.force and os.path.isfile(tgt):
        sContent = stripComments(src)
        tContent = stripComments(tgt)
        if sContent == tContent: return None
    return tgt

def copyFiles(items:set, args:ArgumentParser) -> None:
    cmd = []
    if not args.user: cmd.append(args.sudo)
    cmd.append(args.cp)
    for item in items:
        a = list(cmd)
        a.extend(item)
        logging.info("Copying %s", " ".join(a))
        if not args.dryrun:
            subprocess.run(a, shell=False, check=True)

def mkSystemctl(args:ArgumentParser, options:list=None, extras:set=None, chk:bool=True) -> list:
    cmd = [args.systemctl, "--user"] if args.user else [args.sudo, args.systemctl]
    if options: cmd.extend(options)
    if extras: cmd.extend(extras)
    logging.info("%s", " ".join(cmd))
    if not args.dryrun:
        subprocess.run(cmd, shell=False, check=chk)

def common(args:ArgumentParser) -> tuple:
    if not args.serviceDirectory:
        args.serviceDirectory = "~/.config/systemd/user" if args.user else "/etc/systemd/system"
    args.serviceDirectory = os.path.abspath(os.path.expanduser(args.serviceDirectory))

    services = set()
    timers= set()
    toStart = set()

    for service in args.service:
        service = os.path.abspath(os.path.expanduser(service))
        if not os.path.isfile(service):
            logging.error("%s does not exist", service)
            return 1
        services.add(service)
        dirname = os.path.dirname(service)
        (basename, suffix) = os.path.splitext(os.path.basename(service))
        timer = os.path.join(dirname, basename + ".timer") # Potential timer file
        if os.path.isfile(timer):
            timers.add(timer)
            toStart.add(os.path.basename(timer))
        else:
            toStart.add(os.path.basename(service))

    todos = set(map(os.path.basename, services.union(timers))) # All services and timers basename
    return (services, timers, toStart, todos)

def install(args:ArgumentParser) -> int:
    (services, timers, toStart, todos) = common(args)

    if args.logdir: args.logdir = makeDirectory(args.logdir, args, True) 
    args.serviceDirectory = makeDirectory(args.serviceDirectory, args)

    toCopy = set() # Files that need to be copied
    for fn in services.union(timers): # Copy services and timers as needed
        tgt = needsToBeCopied(fn, args)
        if tgt: toCopy.add((fn, tgt))

    if not toCopy:
        logging.info("Nothing needs to be done")
        return 0

    mkSystemctl(args, ("stop",), toStart, False) # Stop all the processes that need stopped
    mkSystemctl(args, ("disable",), toStart, False) # disable all the services/timers that need stopped
    copyFiles(toCopy, args)
    mkSystemctl(args, ("daemon-reload",)) # Force reload of the daemon
    mkSystemctl(args, ("enable",), todos)

    if toStart: mkSystemctl(args, ("start",), toStart)

    if args.user:
        cmd = (args.loginctl, "enable-linger")
        logging.info("Enable Linger %s", " ".join(cmd))
        if not args.dryrun: subprocess.run(cmd, shell=False, check=True)

    mkSystemctl(args, ("--no-pager", "status"), todos, False)

    if timers:
        mkSystemctl(args, ("--no-pager", "list-timers"), map(os.path.basename, timers), False)

    return 0

def uninstall(args:ArgumentParser) -> int:
    (services, timers, toStart, todos) = common(args)

    toDelete = set()
    for fn in todos:
        ofn = os.path.join(args.serviceDirectory, os.path.basename(fn))
        if os.path.isfile(ofn): toDelete.add(ofn)

    if not toDelete: return 0 # Nothing to be removed
    mkSystemctl(args, ("stop",), toStart, False)
    mkSystemctl(args, ("disable",), todos, False)

    cmd = [args.rm, "-f"] if args.user else [args.sudo, args.rm, "-f"]
    cmd.extend(toDelete)
    subprocess.run(cmd, shell=False, check=True)

    mkSystemctl(args, ("daemon-reload",))
    return 0

def addArgs(parser:ArgumentParser) -> None:
    grp = parser.add_mutually_exclusive_group()
    grp.add_argument("--install", action="store_true", help="Install services and timers")
    grp.add_argument("--uninstall", action="store_true", help="remove services and timers")

    grp = parser.add_mutually_exclusive_group()
    grp.add_argument("--user", action="store_true", help="Install in user space")
    grp.add_argument("--system", action="store_true", help="Install in system space")

    grp = parser.add_argument_group("Command paths")
    grp.add_argument("--sudo", type=str, default="/usr/bin/sudo", help="sudo executable")
    grp.add_argument("--systemctl", type=str, default="/usr/bin/systemctl",
                     help="systemctl executable")
    grp.add_argument("--loginctl", type=str, default="/usr/bin/loginctl",
                     help="loginctl executable")
    grp.add_argument("--mkdir", type=str, default="/bin/mkdir", help="mkdir executable")
    grp.add_argument("--cp", type=str, default="/bin/cp", help="cp executable")
    grp.add_argument("--rm", type=str, default="/bin/rm", help="rm executable")

    grp = parser.add_argument_group("Service/timer related options")
    grp.add_argument("--force", action="store_true", help="Force reloading ...")
    grp.add_argument("--serviceDirectory", type=str, help="Where to copy service file to")
    grp.add_argument("--service", type=str, required=True, action="append", help="Service file(s)")
    grp.add_argument("--logdir", type=str, default="~/logs", help="Where logfiles are stored")
    grp.add_argument("--dryrun", action="store_true", help="Do not actually install anything")

if __name__ == "__main__":
    import Logger

    parser = ArgumentParser()
    Logger.addArgs(parser)
    addArgs(parser)
    args = parser.parse_args()

    Logger.mkLogger(args, fmt="%(asctime)s %(levelname)s: %(message)s")

    sys.exit(uninstall(args) if args.uninstall else install(args))
