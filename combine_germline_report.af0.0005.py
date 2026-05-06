#!/usr/bin/env python3

import os
import glob
import argparse

def main():
    parser = argparse.ArgumentParser(description="Combine TSV files (use header from first file)")
    parser.add_argument("-i", "--input_dir", required=True, help="Input directory")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file")
    args = parser.parse_args()

    files = sorted(glob.glob(os.path.join(args.input_dir, "**", "*iltered.af0.0005.withrc.filtered.tsv"), recursive=True))

    if not files:
        print("No files found!")
        return

    print(f"Found {len(files)} files")

    with open(args.output, "w") as out:
        # --- header from first file ---
        with open(files[0], "r") as f:
            out.write(f.readline())

        # --- append all rows ---
        for f in files:
            if os.path.getsize(f) == 0:
                continue
            with open(f, "r") as infile:
                next(infile)  # skip header
                for line in infile:
                    out.write(line)

    print(f"Done: {args.output}")

if __name__ == "__main__":
    main()
