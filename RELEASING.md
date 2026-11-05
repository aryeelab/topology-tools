# Releasing topology-tools components

This repository holds several independently released components, so release tags are prefixed with
the component name: `<component>/vX.Y.Z` (for example `microc-qc/v1.0.0`). Older repo-wide tags
(`v1.0`, `v1.1`, `v1.1bN`) refer to the WDL pipelines.

Versions follow semantic versioning:
- **patch**: bug fix, no change to the CLI or output schema/values;
- **minor**: new options or output fields, backward compatible; also any change to output values (say so in the notes);
- **major**: breaking change to the CLI or the output schema.

Released versions are referred to by number, never by commit. Tags are annotated and never moved or
deleted; a fix is the next patch release.

## microc-qc (`microc-qc.py`)
1. Tests pass: `python -m pytest tests/test_microc_qc.py` (in the release env).
2. For changes that affect `qc.json`, validate on a real pairs file and compare with the previous release
   (which fields changed and why; peak memory and runtime with `/usr/bin/time -v`).
3. In one release commit:
   - set `__version__` in `microc-qc.py` (it is written to `qc.json` as `microc_qc_version` and printed by `--version`);
   - add the exact environment the release was tested in:
     `envs/microc-qc-<version>.linux-64.explicit.txt` (`conda list -p <env> --explicit`) and
     `envs/microc-qc-<version>.linux-64.pip.txt` (`pip freeze` without `@ file://` lines).
     Its package URLs must start with `https://conda.anaconda.org/`: a conda that installs through a mirror
     (`channel_alias`) writes the mirror's URLs, so replace that prefix before committing.
4. Tag and push:
   ```bash
   git tag -a microc-qc/v<version> -m "microc-qc <version>: <one-line summary>"
   git push origin main microc-qc/v<version>
   ```
   CI (`microc-qc-version-check`) fails if the tag does not match `__version__`.
5. Create the GitHub release from the tag, with notes listing changes since the previous microc-qc tag.

### Installing a released version
```bash
V=1.0.0; SRC=<prefix>/topology-tools/microc-qc-$V; ENV=<envs>/microc-qc-$V
git clone -q https://github.com/aryeelab/topology-tools.git "$SRC" && (cd "$SRC" && git checkout -q microc-qc/v$V)
conda create -q -y -p "$ENV" --file "$SRC/envs/microc-qc-$V.linux-64.explicit.txt"
"$ENV/bin/pip" install -q --no-deps -r "$SRC/envs/microc-qc-$V.linux-64.pip.txt"
(cd "$SRC" && "$ENV/bin/python" -m pytest -q tests/test_microc_qc.py)
"$ENV/bin/python" "$SRC/microc-qc.py" --version
```

## microc-preprocess (`nf/microc_preprocess.nf`, `microc-nf-submit.sh`)
FASTQ → pairs + BAM (see `MICROC-NF.md`). Any change to parse/dedup options, tool versions or MAPQ changes
output values, so it is at least a minor release.
1. Smoke test on the cluster: a few million read pairs with a small `CHUNK_READS` (several chunks, so merge is
   exercised), `OUTDIR` and `DATASET_DIR` in scratch. Check the provenance JSON, the manifest
   row, `tool_versions.txt`, the MAPQ in the pairs header, and compare the pairs by readID with the previous
   release (or with production pairs for the same reads).
2. In one release commit:
   - set `manifest.version` in `nf/nextflow.config` (written to the provenance JSON as `pipeline_version`)
     and `PIPELINE_VERSION` in `microc-nf-submit.sh`;
   - add the exact environment the release was tested in:
     `envs/microc-preprocess-<version>.linux-64.explicit.txt` (`conda list -p <env> --explicit`);
     if a tool version changed, update the version check in `microc-nf-submit.sh`.
     Its package URLs must start with `https://conda.anaconda.org/`: a conda that installs through a mirror
     (`channel_alias`) writes the mirror's URLs, so replace that prefix before committing.
   - update `skills/microc-preprocessing/SKILL.md` (release number, and any change to versions, steps or options).
3. Tag and push:
   ```bash
   git tag -a microc-preprocess/v<version> -m "microc-preprocess <version>: <one-line summary>"
   git push origin main microc-preprocess/v<version>
   ```
4. Check the `microc-preprocess-version-check` CI run (tag = `manifest.version` = `PIPELINE_VERSION`, env file present).
5. Create the GitHub release from the tag (`--latest=false`), with notes listing changes since the previous microc-preprocess tag.

### Installing a released version
```bash
V=1.1.2; SRC=<prefix>/topology-tools/microc-preprocess-$V; ENV=<envs>/microc-preprocess-$V
git clone -q https://github.com/aryeelab/topology-tools.git "$SRC" && (cd "$SRC" && git checkout -q microc-preprocess/v$V)
conda create -q -y -p "$ENV" --file "$SRC/envs/microc-preprocess-$V.linux-64.explicit.txt"
"$ENV/bin/pairtools" --version && JAVA_HOME="$ENV" "$ENV/bin/nextflow" -version
```
Run `bash $SRC/microc-nf-submit.sh` from anywhere. Nextflow runs in `OUTDIR` (its cache `.nextflow/`,
`.nextflow.log`, the master job log `<SAMPLEID>.master.<jobid>.log` and the per-launch trace `<SAMPLEID>.trace.<timestamp>.txt` are there; the checkout is never written to)
and its work dir is `OUTDIR/work`; pairs and BAM are copied out, so `work/` can be removed once a
run's outputs are in place.

## preprocess-hic (`nf/preprocess-hic/`, `preprocess-hic-submit.sh`, `preprocess-hic-pool.sh`)
Hi-C FASTQ → lossless BAM, pairs, QC, mcool (see `PREPROCESS-HIC.md`). Any change to align/parse/dedup/restrict/cooler
options or tool versions changes output values, so it is at least a minor release.
1. Tests pass in the release env, with its `bin` on `PATH` (the end-to-end test runs Nextflow locally on synthetic data):
   `PATH=<env>/bin:$PATH JAVA_HOME=<env> python -m pytest -v tests/preprocess-hic`
2. Smoke test on the cluster: a few million read pairs of a real library with a small `CHUNK_READS` (several chunks, so
   merge is exercised), `OUTDIR` in scratch. Check the provenance JSON, `hic_qc.json` (lossless check `pass`, mapping,
   duplication, cis/trans and fragment-artefact rates plausible for the protocol), the parse options in the pairs header,
   and compare with the previous release's outputs for the same reads (pairs by readID, stats, mcool sums).
3. In one release commit:
   - set `manifest.version` in `nf/preprocess-hic/nextflow.config` (written to the provenance JSON, `hic_qc.json` and the
     cooler metadata) and `PIPELINE_VERSION` in `preprocess-hic-submit.sh` and `preprocess-hic-pool.sh`;
   - add the exact environment the release was tested in: `envs/preprocess-hic-<version>.linux-64.explicit.txt`
     (`conda list -p <env> --explicit`), with `https://conda.anaconda.org/` URLs (replace a mirror's prefix);
     if a tool version changed, update the checks in `preprocess-hic-submit.sh`.
   Put the smoke-test numbers in the commit message.
4. Tag and push:
   ```bash
   git tag -a preprocess-hic/v<version> -m "preprocess-hic <version>: <one-line summary>"
   git push origin main preprocess-hic/v<version>
   ```
5. Check the `preprocess-hic-version-check` CI run; create the GitHub release from the tag (`--latest=false`).

### Installing a released version
```bash
V=1.0.0; SRC=<prefix>/topology-tools/preprocess-hic-$V; ENV=<envs>/preprocess-hic-$V
git clone -q https://github.com/aryeelab/topology-tools.git "$SRC" && (cd "$SRC" && git checkout -q preprocess-hic/v$V)
conda create -q -y -p "$ENV" --file "$SRC/envs/preprocess-hic-$V.linux-64.explicit.txt"
(cd "$SRC" && PATH="$ENV/bin:$PATH" JAVA_HOME="$ENV" "$ENV/bin/python" -m pytest -q tests/preprocess-hic)
```
