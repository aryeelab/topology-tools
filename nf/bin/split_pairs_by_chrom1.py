#!/usr/bin/env python3
"""Split a sorted .pairs stream (stdin) into one bgzipped file per chrom1, each with the full header.

Usage: split_pairs_by_chrom1.py PREFIX [BGZIP_THREADS]

Writes PREFIX.blk_<id>.pairs.gz for each chrom1 and PREFIX.blocks.tsv (<id> TAB <chrom1>).
<id> is the chromosome name when it only has [A-Za-z0-9._-], otherwise "x" + its hex encoding
(e.g. "!" for unmapped sides becomes "x21"). The input must be sorted by chrom1 (pairtools sort),
so each chrom1 is one contiguous block; a chrom1 that comes back later is an error.
"""
import re
import subprocess
import sys

SAFE = re.compile(rb"^[A-Za-z0-9._-]+$")


def block_id(chrom: bytes) -> str:
    return chrom.decode() if SAFE.match(chrom) else "x" + chrom.hex()


def main() -> None:
    prefix = sys.argv[1]
    threads = sys.argv[2] if len(sys.argv) > 2 else "2"
    header = []
    seen = []
    chrom, proc, buf = None, None, []

    def close_block():
        if proc is not None:
            proc.stdin.writelines(buf)
            proc.stdin.close()
            if proc.wait() != 0:
                sys.exit(f"bgzip failed for block {chrom!r}")

    for line in sys.stdin.buffer:
        if line.startswith(b"#"):
            if proc is not None:
                sys.exit("header line after data lines")
            header.append(line)
            continue
        c = line.split(b"\t", 2)[1]
        if c != chrom:
            close_block()
            if c in seen:
                sys.exit(f"input not sorted by chrom1: {c!r} seen twice")
            chrom = c
            seen.append(c)
            proc = subprocess.Popen(["bgzip", "-c", "-@", threads],
                                    stdin=subprocess.PIPE,
                                    stdout=open(f"{prefix}.blk_{block_id(c)}.pairs.gz", "wb"))
            proc.stdin.writelines(header)
            buf = []
        buf.append(line)
        if len(buf) >= 100000:
            proc.stdin.writelines(buf)
            buf = []
    close_block()
    with open(f"{prefix}.blocks.tsv", "w") as f:
        for c in seen:
            f.write(f"{block_id(c)}\t{c.decode()}\n")


if __name__ == "__main__":
    main()
