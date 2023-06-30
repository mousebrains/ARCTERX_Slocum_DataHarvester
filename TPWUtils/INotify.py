#! /usr/bin/env python3
#
# Use pyinotify to handle file notification and send
# the notification to a queue
#

from argparse import ArgumentParser
import pyinotify
import queue
import logging
import time
import os
try:
    from .Thread import Thread # As a module
except:
    from Thread import Thread # From within module

class INotify(Thread):
    def __init__(self, args:ArgumentParser, flags:int=None) -> None:
        Thread.__init__(self, "INotify", args)
        self.__wm = pyinotify.WatchManager()
        self.__notifier = pyinotify.Notifier(self.__wm)
        self.queue = queue.Queue()
        if flags is not None:
            self.__flags = flags
        else:
            self.__flags = pyinotify.IN_CREATE 
            self.__flags|= pyinotify.IN_MODIFY 
            self.__flags|= pyinotify.IN_CLOSE_WRITE 
            self.__flags|= pyinotify.IN_MOVED_TO 
            self.__flags|= pyinotify.IN_MOVED_FROM
            self.__flags|= pyinotify.IN_MOVE_SELF
            self.__flags|= pyinotify.IN_DELETE
            self.__flags|= pyinotify.IN_DELETE_SELF

    @staticmethod
    def __maskname(mask:int) -> str:
        # This should be able to use pyinotify.EventsCodes.maskname, but it fails
        items = []
        codes = pyinotify.EventsCodes.FLAG_COLLECTIONS["OP_FLAGS"]
        for key in codes:
            if mask & codes[key]: items.append(key)
        return "|".join(items) if items else None

    def addTree(self, tgt:str) -> None:
        self.addWatch(tgt, qRecursive=True, qAutoAdd=True)

    def addWatch(self, tgt:str, mask:int=None, qRecursive:bool=False, qAutoAdd:bool=False) -> bool:
        tgt = os.path.abspath(os.path.expanduser(tgt))
        if os.path.isdir(tgt):
            mask = mask if mask is not None else self.__flags
            self.__wm.add_watch(path=tgt, mask=mask, proc_fun=self.__eventHandler, 
                    rec=qRecursive, auto_add=qAutoAdd)
            logging.info("Added watch for %s, rec %s auto %s msk %s",
                    tgt, qRecursive, qAutoAdd, self.__maskname(mask))
            return True
        logging.error("Path %s does not exist", tgt)
        return False

    def runIt(self) -> None: # Called on thread start
        logging.info("Starting loop")
        self.__notifier.loop() # All the action happens in __eventHandler
        logging.warning("Leaving loop")

    def __eventHandler(self, e:pyinotify.Event) -> None:
        t0 = time.time() # Time of the event
        fn = e.path if e.dir else os.path.join(e.path, e.name)
        self.queue.put((t0, fn))
        logging.debug("Event %s, %s", fn, e.maskname)

if __name__ == "__main__":
    import Logger

    class Reader(Thread):
        def __init__(self, args:ArgumentParser, q:queue.Queue) -> None:
            Thread.__init__(self, "Reader", args)
            self.__queue = q

        def runIt(self) -> None:
            q = self.__queue
            while True:
                (t0, fn) = q.get()
                q.task_done()
                logging.info("%s %s", t0, fn)

    parser = ArgumentParser()
    Logger.addArgs(parser)
    parser.add_argument("tgt", nargs="+", help="Directories to watch")
    args = parser.parse_args()

    Logger.mkLogger(args)

    i = INotify(args)
    rdr = Reader(args, i.queue)
    i.start()
    rdr.start()
    for tgt in args.tgt: i.addTree(tgt)

    try:
        Thread.waitForException()
    except:
        logging.exception("Exception from INotify")
