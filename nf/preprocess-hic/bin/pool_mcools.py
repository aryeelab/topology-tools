#!/usr/bin/env python3
"""Pool replicate .mcool files from preprocess-hic into one balanced .mcool.

Usage: pool_mcools.py --pool-id ID --outdir DIR [--nproc N] [--git-commit C] [--version V] A.mcool B.mcool [...]

Checks that the inputs were made by preprocess-hic with the same MAPQ filter (cooler metadata) and have the
same bins and resolutions, sums them at the base (smallest) resolution with `cooler merge`, then coarsens and
balances with `cooler zoomify --balance` at the inputs' resolutions. Writes DIR/ID.mapq_<Q>.mcool (never
overwrites) and DIR/ID.mapq_<Q>.pool.provenance.json.
"""
import argparse
import datetime
import json
import os
import socket
import subprocess
import sys

import cooler
import h5py


def run(cmd, log):
    log.append(cmd)
    print("+", " ".join(cmd), file=sys.stderr, flush=True)
    subprocess.run(cmd, check=True)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--pool-id", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--nproc", type=int, default=4)
    ap.add_argument("--git-commit", default=None)
    ap.add_argument("--version", default=None)
    ap.add_argument("mcools", nargs="+")
    a = ap.parse_args(argv)
    started = datetime.datetime.now().astimezone().isoformat(timespec="seconds")
    if len(a.mcools) < 2:
        sys.exit("pool_mcools.py: give at least two .mcool files")

    inputs, ress, metas, chroms = [], [], [], []
    for m in a.mcools:
        m = os.path.abspath(m)
        res = sorted(int(p.rsplit("/", 1)[1]) for p in cooler.fileops.list_coolers(m))
        c = cooler.Cooler(f"{m}::/resolutions/{res[0]}")
        meta = c.info.get("metadata") or {}
        if meta.get("pipeline") != "preprocess-hic" or "min_mapq" not in meta:
            sys.exit(f"pool_mcools.py: {m} has no preprocess-hic metadata (pipeline, min_mapq)")
        st = os.stat(m)
        inputs.append({"path": m, "size": st.st_size,
                       "mtime": datetime.datetime.fromtimestamp(st.st_mtime).astimezone().isoformat(timespec="seconds"),
                       "sample_id": meta.get("sample_id"), "min_mapq": meta["min_mapq"],
                       "pipeline_version": meta.get("pipeline_version"), "sum_base": int(c.info["sum"])})
        ress.append(res)
        metas.append(meta)
        chroms.append(list(zip(c.chromnames, map(int, c.chromsizes))))
    if len({i["min_mapq"] for i in inputs}) != 1:
        sys.exit(f"pool_mcools.py: inputs differ in min_mapq: {[i['min_mapq'] for i in inputs]}")
    if any(r != ress[0] for r in ress):
        sys.exit(f"pool_mcools.py: inputs differ in resolutions: {ress}")
    if any(c != chroms[0] for c in chroms):
        sys.exit("pool_mcools.py: inputs differ in chromosomes or sizes")

    q, base = inputs[0]["min_mapq"], ress[0][0]
    os.makedirs(a.outdir, exist_ok=True)
    out = os.path.join(a.outdir, f"{a.pool_id}.mapq_{q}.mcool")
    if os.path.exists(out):
        sys.exit(f"pool_mcools.py: {out} exists; not overwriting")
    tmp = os.path.join(a.outdir, f".{a.pool_id}.mapq_{q}.base.cool")
    log = []
    run(["cooler", "merge", tmp] + [f"{i['path']}::/resolutions/{base}" for i in inputs], log)
    meta = {"pipeline": "preprocess-hic", "pipeline_version": a.version, "pool_id": a.pool_id, "min_mapq": q,
            "pooled_from": [i["sample_id"] for i in inputs]}
    with h5py.File(tmp, "r+") as f:
        f.attrs["metadata"] = json.dumps(meta)
    run(["cooler", "zoomify", "-p", str(a.nproc), "--balance", "-r", ",".join(map(str, ress[0])), "-o", out + ".part", tmp], log)
    os.rename(out + ".part", out)
    os.remove(tmp)
    pooled_sum = int(cooler.Cooler(f"{out}::/resolutions/{base}").info["sum"])
    expected = sum(i["sum_base"] for i in inputs)
    prov = {
        "pipeline": "preprocess-hic (pool_mcools.py)", "pipeline_version": a.version, "git_commit": a.git_commit,
        "pool_id": a.pool_id, "min_mapq": q, "resolutions": ress[0], "inputs": inputs, "output": out,
        "sum_base": pooled_sum, "sum_check": "pass" if pooled_sum == expected else f"FAIL (expected {expected})",
        "commands": [" ".join(c) for c in log],
        "tool_versions": {"cooler": cooler.__version__, "python": sys.version.split()[0]},
        "slurm_job": os.environ.get("SLURM_JOB_ID"), "node": socket.gethostname(), "started": started,
        "completed": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
    }
    with open(os.path.join(a.outdir, f"{a.pool_id}.mapq_{q}.pool.provenance.json"), "w") as f:
        f.write(json.dumps(prov, indent=2) + "\n")
    if pooled_sum != expected:
        sys.exit(f"pool_mcools.py: pooled sum {pooled_sum} != sum of inputs {expected}")


if __name__ == "__main__":
    main()
