#!/usr/bin/env python3
# Author: Zifo Bioinformatics
# Email: bioinformatics@zifornd.com
# License: MIT

import itertools
import yaml
import platform


def main(input_file: str, output_file: str) -> None:

    # Open file
    with open(input_file, "r") as h:
        header = h.readline()

    # Split header by delimiter (e.g., transcript_GBR_1)
    samples = header.split("\\t")

    # Trim replicate number (e.g., transcript_GBR_1 -> transcript_GBR)
    conditions = [sample.rsplit("_", 1)[0] for sample in samples]

    # Create list of consecutive condition indices (e.g., [[1,2], [3,4]])
    last_index = 0
    out = []
    for v, g in itertools.groupby(enumerate(conditions), lambda k: k[1]):
        l = [*g]
        out.append([last_index + 1, l[-1][0] + 1])
        last_index += len(l)

    # Assert that condition indices are consecutive
    if len(out) != 2:
        raise ValueError("Column numbers have to be continuous, with no overlapping or missing columns between them.")

    # Format ranges for printing
    groups = ",".join([f"{start}-{end}" for start, end in out])

    # Write to output file
    with open(output_file, "w") as h:
        h.write(f"{groups}\\n")

if __name__ == "__main__":
    prefix = "${task.ext.prefix}" if "${task.ext.prefix}" not in ["null", ""] else "${cond1}-${cond2}"
    input_file = "${psivec}"
    output_file = f"{prefix}_groups.txt"
    main(input_file, output_file)

    # Write version information to versions.yml
    versions = {"${task.process}": {"python": platform.python_version()}}

    with open("versions.yml", "w") as f:
        yaml.dump(versions, f)
