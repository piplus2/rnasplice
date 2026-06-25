#!/usr/bin/env python
# Author: Zifo Bioinformatics
# Email: bioinformatics@zifornd.com
# License: MIT

## Update - author: Paolo Inglese [2026-06-23]
# Fix regex expression and improved the file management

import argparse
import glob
import os
import re
import shelve


def main(prefix):
    # Check for non alphanumerics
    regex = re.compile("[@!#$%^&*()<>?/\\|}{~:]")
    if regex.search(prefix) is not None:
        raise ValueError(
            "Prefix contains special chars: " + prefix + ". Please remove special chars: [@!#$%^&*()<>?/\\|}{~:]"
        )

    # Open and extract shelve file contents
    shelve_file = shelve.open(prefix + "/genes_to_filenames.shelve")
    data = [(gene, filename) for gene, filename in shelve_file.items()]
    shelve_file.close()

    # Remove existing shelve files
    for filename in glob.glob(prefix + "/genes_to_filenames.shelve.*"):
        os.remove(filename)

    # Create new shelve file with relative paths
    shelve_file = shelve.open(prefix + "/genes_to_filenames.shelve")
    for gene, filename in data:
        elements = filename.split(prefix)
        print(elements)
        if len(elements) != 2:
            raise ValueError(
                "Filename contains " + str(len(elements) - 1) + " instances of index delimiter. Expected only one."
            )
        shelve_file[gene] = "./" + prefix + elements[1]
    shelve_file.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse Miso index shelve files")
    parser.add_argument(
        "-p",
        "--prefix",
        default="index",
        type=str,
        help="Folder where index shelve files stored",
    )

    args = parser.parse_args()
    main(args.prefix)
