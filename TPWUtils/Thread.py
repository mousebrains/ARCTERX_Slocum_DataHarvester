#! /usr/bin/env python3
#
# A thread super class that handles exceptions from the actual thread.
# The actual thread run method is called runIt, otherwise it is like a normal
# thread class.
#
# June-2021, Pat Welch, pat@mousebrains.com

from argparse import ArgumentParser
import threading
import queue
#
# Base class for threading which catches exceptions and sends them to a queue
#
class Thread(threading.Thread):
    '''
    A thread super class that handles exceptions from the actual thread.
    The actual thread run method is called runIt.
    Otherwise the actual thread behaves like a "normal" thread
    '''

    __queue = queue.Queue() # Static variable

    def __init__(self, name:str, args:ArgumentParser=None) -> None:
        '''
        name: is the name of the thread saved in self.name and used by logging messages
        args: is saved in self.args
        '''
        threading.Thread.__init__(self, daemon=True)
        self.name = name
        self.args = args

    def run(self) -> None: # Called on thread start
        try:
            self.runIt() # Call the actual class's run function inside a try stanza
        except Exception as e:
            self.__queue.put(e)

    @classmethod
    def isQueueEmpty(cls) -> bool: 
        return cls.__queue.empty()

    @classmethod
    def waitForException(cls, timeout=None):
        if timeout is None:
            e = cls.__queue.get()
            raise e
        try:
            e = cls.__queue.get(timeout=timeout)
            print(e)
        except queue.Empty:
            return
        except Exception as e:
            raise e

if __name__ == "__main__":
    import Logger
    import logging
    import time

    class A(Thread):
        def __init__(self, args:ArgumentParser):
            Thread.__init__(self, "A", args)

        @staticmethod
        def addArgs(paresr:ArgumentParser):
            parser.add_argument("--dt", type=float, default=1.7,
                    help="Time to wait to throw an exception")
        def runIt(self):
            ''' I'll throw an exception after --dt seconds '''
            dt = self.args.dt
            logging.info("Going to throw an error after %s seconds", dt)
            time.sleep(dt)
            logging.warning("Throwing a NotImplemented exception")
            raise NotImplemented

    parser = ArgumentParser()
    Logger.addArgs(parser)
    A.addArgs(parser)
    args = parser.parse_args()

    Logger.mkLogger(args)

    thrd = A(args)
    thrd.start()
    try:
        Thread.waitForException()
    except:
        logging.exception("Exception from A")
