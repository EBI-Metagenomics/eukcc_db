#!/usr/bin/env python3

from ftplib import FTP
import ftplib
import os
from urllib.request import urlopen
import hashlib
from functools import partial
import urllib.request
import logging
import argparse

logging.basicConfig(format="%(asctime)s %(message)s", datefmt="%m/%d/%Y %I:%M:%S %p", level=logging.DEBUG)


# reusable class for fetching GCA
class downloadGCA:
    def __init__(self, GCA, directory=".", prot=False):
        # first avlidate input
        # validate GCA
        self.validateGCA(GCA)
        self.prot = prot

        localPath = self.defineLocalPath(directory, GCA)

        # download only is not there
        if not self.checkLocally(localPath):
            # open ftp connection
            ftp = FTP("ftp.ncbi.nlm.nih.gov")
            ftp.login()
            # construct url and fetch file. The control connection is only needed for the
            # listing: NCBI drops it after 60 s idle, i.e. during any longer download
            path = self.constructURL(ftp, GCA)
            ftp.quit()
            # the file itself over HTTPS: urllib's ftp:// handler returned a corrupted
            # 319 MB genome (GCA_046563775.1) without an error; HTTPS checks the length
            url = "https://ftp.ncbi.nlm.nih.gov/{}".format(path)
            self.download(localPath, url, GCA)
        else:
            logging.info("{}: file already downloaded".format(GCA))

    def validateGCA(self, GCA):
        if len(GCA) == 0:
            raise Exception("GCA is empty")
        return True

    def defineLocalPath(self, directory, GCA):
        if not self.prot:
            path = "{}/{}.fna.gz".format(directory, GCA)
        else:
            path = "{}/{}_prot.faa.gz".format(directory, GCA)
        return path

    def checkLocally(self, path):
        return os.path.isfile(path)

    def constructURL(self, ftp, GCA):

        # ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/987/745/GCA_000987745.1_ASM98774v1/
        base = "ftp://ftp.ncbi.nlm.nih.gov/"

        # url = '{}/GCA/{}/{}/{}/'.format(base, GCA[4:7],GCA[7:10], GCA[10:13])
        ftp_path = "{}/{}/{}/{}/{}/".format("genomes/all", GCA[0:3], GCA[4:7], GCA[7:10], GCA[10:13])
        # get directory with full name of GCA
        files = self.listFTP(ftp, ftp_path)
        if len(files) == 1:
            path = files[0]
        else:
            # in case of more than one assembly we search for the
            # GCA
            s = "{}_".format(GCA)
            found = False
            print(files, s)
            for file in files:
                if file.find(s) != -1:
                    path = file
                    found = True
            if not found:
                raise Exception("Could not find any genomic assembly for {}. Although there where many".format(GCA))

        # else go one deeper
        files = self.listFTP(ftp, path)
        # find assembly
        found = False
        # get dirname for fetching the assembly
        dirname = os.path.basename(path)
        if not self.prot:
            # we want to fetch the genome
            fname = "{}_genomic.fna.gz".format(dirname)
        else:
            # we want to fetch the proteome
            fname = "{}_protein.faa.gz".format(dirname)

        for file in files:
            if file.find(fname) != -1:
                path = file
                found = True
                break
        if not found:
            logging.warning("Looked online in folder {}{}".format(base, path))
            raise Exception("Could not find any genomic assembly for {}".format(GCA))

        # now construct url
        return path

    def listFTP(self, ftp, path):
        try:
            files = ftp.nlst(path)
        except ftplib.error_perm:
            raise
        return files

    def download(self, local, remote, GCA):

        # handle = open(local, 'wb')
        # ftp.retrbinary('RETR %s' % remote, handle.write)

        logging.info("Downloading from {}".format(remote))
        urllib.request.urlretrieve(remote, local)
        logging.info("Download done")

        # get md5 checks; a failed check is an error, so Nextflow retries the download
        if self.verifyFile(local, remote, GCA):
            logging.info("{}: fetched assembly".format(GCA))
        else:
            os.remove(local)
            raise Exception("{}: md5 mismatch, file removed".format(GCA))
        return True

    def verifyFile(self, local, remote, GCA):
        logging.info("Fetching md5sum info")
        # md5checksums.txt sits next to the file: no FTP listing needed
        url = "{}/md5checksums.txt".format(os.path.dirname(remote))
        data = urlopen(url)  # it's a file like object and works just like a file
        checksum = "unknown"
        # fetch checksums
        for line in data:  # files are iterable
            s = line.split()
            if s[1].decode("utf-8").find(os.path.basename(remote)) != -1:
                checksum = s[0].decode("utf-8")
                # print(s)
        # compute file checksum
        md5result = self.md5sum(local)
        if md5result != checksum:
            logging.warning("{}: expected {} found {}".format(GCA, checksum, md5result))
        return md5result == checksum

    def md5sum(self, filename):
        d = hashlib.md5()
        with open(filename, "rb") as f:
            for buf in iter(partial(f.read, 1 << 20), b""):
                d.update(buf)
        return d.hexdigest()


def main():
    parser = argparse.ArgumentParser(description="Download content from NCBI")
    parser.add_argument("identifier", type=str, help="Pass GCA or GCF accession")
    parser.add_argument("outdir", type=str, help="specify outdir")
    parser.add_argument("--prot", action="store_true", default=False, help="Download the proteome not the genome")
    args = parser.parse_args()

    downloadGCA(args.identifier, args.outdir, args.prot)


if __name__ == "__main__":
    main()
