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
