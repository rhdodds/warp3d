#!/usr/bin/env python3

import sys

def main():
    if len(sys.argv) != 2:
        print("usage: sum_walltimes.py transcript_file")
        sys.exit(1)

    total = 0.0
    count = 0

    with open(sys.argv[1], "r") as f:
        for line in f:
            if "... Ran:" in line and "walltime:" in line:
                fields = line.split()
                jobname = fields[2]
                wtime = float(fields[-1])
                count += 1
                total += wtime
                print(f"{count:3d}. {jobname:12s}  {wtime:8.2f} sec")

    print(f"\nTotal wall time = {total:.2f} sec")
    print(f"                = {total/60.0:.2f} min")
    print(f"                = {total/3600.0:.3f} hr")

if __name__ == "__main__":
    main()
    