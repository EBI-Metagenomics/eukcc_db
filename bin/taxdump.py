#!/usr/bin/env python3
"""The NCBI taxonomy dump in either archive the FTP serves: taxdump.tar.gz (the dump of
the day, pub/taxonomy/) or taxdmp_YYYY-MM-01.zip (monthly snapshots, taxdump_archive/,
zip only). ete3 only reads a tar.gz.

    taxdump.py date <dump>              date of nodes.dmp in the dump (YYYY-MM-DD)
    taxdump.py extract <dump> <member>  writes that member to the working directory
    taxdump.py targz <dump> <out>       the dump as a tar.gz with the .dmp files at its root
"""
import datetime
import os
import shutil
import struct
import sys
import tarfile
import zipfile


def zip_mtime(info):
    """UTC mtime of a zip member: the extended timestamp field (0x5455) when there is one.
    The DOS time is local to whoever zipped it: NCBI's 2022-04-01 03:27 UTC reads 03-31."""
    extra, i = info.extra, 0
    while i + 4 <= len(extra):
        tag, size = struct.unpack("<HH", extra[i:i + 4])
        if tag == 0x5455 and size >= 5 and extra[i + 4] & 1:
            return struct.unpack("<I", extra[i + 5:i + 9])[0]
        i += 4 + size
    return datetime.datetime(*info.date_time).timestamp()


def date(dump):
    if zipfile.is_zipfile(dump):
        mtime = zip_mtime(zipfile.ZipFile(dump).getinfo("nodes.dmp"))
    else:
        with tarfile.open(dump) as t:
            mtime = t.getmember("nodes.dmp").mtime
    return datetime.datetime.utcfromtimestamp(mtime).strftime("%Y-%m-%d")


def extract(dump, member):
    if zipfile.is_zipfile(dump):
        zipfile.ZipFile(dump).extract(member)
    else:
        with tarfile.open(dump) as t:
            t.extract(member)


def targz(dump, out):
    if not zipfile.is_zipfile(dump):
        shutil.copyfile(dump, out)
        return
    z = zipfile.ZipFile(dump)
    with tarfile.open(out, "w:gz") as t:
        for info in z.infolist():
            if not info.filename.endswith(".dmp"):
                continue
            member = tarfile.TarInfo(os.path.basename(info.filename))
            member.size = info.file_size
            member.mtime = zip_mtime(info)
            with z.open(info) as fin:
                t.addfile(member, fin)


if __name__ == "__main__":
    cmd, args = sys.argv[1], sys.argv[2:]
    if cmd == "date":
        print(date(*args))
    elif cmd == "extract":
        extract(*args)
    elif cmd == "targz":
        targz(*args)
    else:
        sys.exit(__doc__)
