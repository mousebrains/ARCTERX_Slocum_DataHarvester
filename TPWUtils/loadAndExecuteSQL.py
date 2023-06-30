#! /usr/bin/env python3
#
# Read a file of SQL commands and execute it
#
# This has been tested with PostgreSQL 14 and psycopg 3, 
# but it should work with most databases
#
# April-2023, Pat Welch, pat@mousebrains.com

import logging

def loadAndExecuteSQL(db, fn:str, tableName:str=None) -> bool:
    body = None

    try:
        cur = db.cursor()

        if tableName is not None: # Check if this table already exists
            cur.execute("SELECT EXISTS(SELECT relname FROM pg_class WHERE relname=%s);",
                        (tableName,))
            for row in cur:
                if row[0]: return True # Already exists
                break

        with open(fn, "r") as fp: body = fp.read()
        logging.info("Loaded %s, %s bytes", fn, len(body))

        cur.execute("BEGIN TRANSACTION;")
        cur.execute(body)
        db.commit()
        return True
    except:
        logging.exception("Unable to execute %s", fn)
        db.rollback()
        return False

if __name__ == "__main__":
    from argparse import ArgumentParser
    import Logger
    import psycopg

    parser = ArgumentParser()
    Logger.addArgs(parser)
    parser.add_argument("db", type=str, help="Database name")
    parser.add_argument("sql", type=str, help="File containing SQL statements")
    args = parser.parse_args()


    Logger.mkLogger(args, fmt="%(asctime)s %(levelname)s: %(message)s")

    with psycopg.connect(f"dbname={args.db}") as db:
        loadAndExecuteSQL(db, args.sql)
