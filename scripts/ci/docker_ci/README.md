# OpenIFS Docker CI test — control vs test branch

Builds two Docker images (one per branch). In each container, runs
`openifs-test.sh -cb` (configure + build), then `openifs-test.sh -t`
(ctest) with `IFS_TEST_BITIDENTICAL=init IFS_TEST_LEGACY=1` so the
framework drops a `SAVED_NORMS` reference file in every `test*/` subdir.
The two stages are run separately so build and ctest output can be captured
in isolation as separate artifacts.

The bit-comparison runs **inside** the test container — the host never
reads NORMS data. The host receives a text report
(`<control>-<sha7>__<test|dir>.txt`) plus raw per-branch build, ctest, and
LastTest output files. Three exit codes: **0** PASS, **1** FAIL, **2**
INCONCLUSIVE (control failed, bit-comparison skipped — see below).

## Prerequisites

- Docker installed and running
- Python 3 with `gitpython` and `pyyaml` (same venv as `../docker/`)
- Git configured with SSH access to the OpenIFS repository

```bash
cd scripts/ci/docker_ci

# Re-use the venv from the sibling docker/ builder
source ../docker/openifs-env/bin/activate     # (or create a new one)
python3 -m pip install gitpython pyyaml
```

## Quick start

1. Edit `config/ci_test_docker.yml` — set `control_branch` and
   `test_branch` (and `openifs_repo_url` if not the default).

2. Run:

   ```bash
   python3 ci-oifs-docker.py -c config/ci_test_docker.yml
   ```

3. Read the report at `<ci_reports>/<control>-<sha7>__<test|dir>.txt` (the filename encodes the control branch + its 7-char commit SHA, and either the test branch or the basename of the local source dir when `test_branch` is empty or a path). The script
   exits **0** on PASS, **1** on FAIL, **2** on INCONCLUSIVE.

## What runs where

```
host                     control container        test container
────────────────────────────────────────────────────────────────────────────
ci-oifs-docker.py
  ├── build & test ─────► openifs-test.sh -cb, then -t (output captured)
  │   ◄── docker cp ──── ctest output  ──► openifs_test_output_control.txt
  │   ◄── docker cp ──── SAVED_NORMS tarball
  │   (control failure here is tolerated → INCONCLUSIVE; script continues)
  │
  ├── build & test ───────────────────────────────► openifs-test.sh -cb, then -t (captured)
  │   ◄── docker cp ─────────────────────────────── ctest output  ──► openifs_test_output_test.txt
  │
  ├── docker cp control tarball + bitcompare.py INTO test container
  ├── docker exec ────────────────────────────────► python3 bitcompare --report
  │   ◄── docker cp ─────────────────────────────── <control>-<sha7>__<test|dir>.txt
  │
   ├── append CI summary to the host-side report
  │
  └── exit 0 (PASS) | 1 (FAIL) | 2 (INCONCLUSIVE — control failed)
```

## Output and exit codes

After a run, these artefacts land in `<ci_reports>` when available:

| File | Contents |
| ---- | -------- |
| `<control>-<sha7>__<test\|dir>.txt` | Bit-comparison or skipped-comparison report plus CI summary |
| `openifs_build_output_control.txt` | Raw configure/build stdout/stderr from the control container |
| `openifs_test_output_control.txt` | Raw ctest stdout/stderr from the control container |
| `openifs_lasttest_output_control.txt` | CTest `LastTest.log` from the control container |
| `openifs_build_output_test.txt` | Raw configure/build stdout/stderr from the test container |
| `openifs_test_output_test.txt` | Raw ctest stdout/stderr from the test container |
| `openifs_lasttest_output_test.txt` | CTest `LastTest.log` from the test container |

The `<control>-<sha7>__…` file is the headline artefact. Its sections, in order:

1. **Bit-comparison body** — per-test PASS/FAIL lines from
   `openifs_branch_bitcompare.py`, ending with `RESULT: PASS|FAIL`.
   Replaced by a `BIT-COMPARISON SKIPPED` header + `RESULT: SKIPPED` when
   the control phase failed and there were no SAVED_NORMS to compare.
2. **`CI SUMMARY`** — branches with annotation
   (`built + tested` / `reused cached NORMS` / `FAILED — bit-comparison skipped`),
   per-branch build and test status, bit-comparison verdict and skip reason,
   final result (`PASS` / `FAIL` / `INCONCLUSIVE`), report path, control
   NORMS path, per-phase timings, total wall-clock.

### Exit codes

| Code | Meaning |
| ---- | ------- |
| 0    | **PASS** — both branches built and ctest'd, all SAVED_NORMS bit-matched |
| 1    | **FAIL** — both branches built and ctest'd, at least one test's SAVED_NORMS differed |
| 2    | **INCONCLUSIVE** — control phase failed (build or ctest); bit-comparison skipped. Useful when the test branch is the fix for a broken main: the test side still builds and runs, you just can't bit-compare against a missing reference. The test branch's own output files are still uploaded as separate artifacts |

The **test phase is strict**: a test-branch build/ctest failure is a hard
error and aborts the script (no INCONCLUSIVE fallback for that side).
For test-side debugging, set `remove_test_container: False` and inspect
the container with `docker exec -it oifs-ci-test /bin/bash`.

## Key config flags

| Flag | Purpose |
| ---- | ------- |
| `control_branch`                          | Remote branch name to use as the reference (always cloned, unless `clone_openifs_control: False` and you pre-populate the build dir) |
| `test_branch`                             | Accepts: empty (auto-resolves to `<ci-oifs-docker.py>/../..`), a directory path (used as a local source — `~` is expanded), or a remote branch name (cloned from `openifs_repo_url`). Detected automatically. Local sources are copied into the build dir, skipping `.git`, `build/`, `.cache`, `.bootstrap`, `__pycache__`, `openifs-env` |
| `openifs_repo_url`                        | Source repo (SSH) |
| `base_docker_image`                       | GCC tag from docker.io/library/gcc |
| `openifs_test_extra_flags`                | Extra configure-time args forwarded to `openifs-test.sh -cb` (e.g. `--without-double-precision --cmake=BUILD_ifsbench=OFF --clean`). Not passed to the `-t` ctest call |
| `clone_openifs_control`                   | When False, do not clone the control branch — expect the source pre-populated at `<openifs_build_docker_dir>/build_dir_control/<openifs_version>/`. Test has no equivalent flag (auto-detected from `test_branch`) |
| `reuse_control_if_present`                | Skip the entire control phase if the image-keyed control NORMS tarball (`control_saved_norms_<version>_gcc<image>.tgz`) already exists |
| `remove_test_container`                   | Whether to `docker rm -f` containers on success |
| `force_rebuild`, `force_reclone`          | Same semantics as `../docker/config/create_openifs_docker.yml`. Recommend `force_rebuild: False` for fast iteration: control and test images share their first ~80% of layers (apt install, OpenMPI build, user setup), so the second sequential build reuses the first's cache and only re-runs from `COPY ${OPENIFS_DIR}` onward. Flip to `True` only if you suspect cache corruption |

## Files

| File                              | Role |
| --------------------------------- | ---- |
| `ci-oifs-docker.py`               | Driver — orchestrates both branches and the in-container comparison |
| `Dockerfile.ci`                   | Minimal image — only the toolchain (apt + OpenMPI build) and the OpenIFS source. No SCM data, no experiment data, no ifsdata/rtables/climate (ctest fetches its own inputs from ECPDS via `ifs-test/bin/storage.py`) |
| `openifs_branch_bitcompare.py`    | SAVED_NORMS comparator (per-test `filecmp`-based, with `--report`) |
| `config/ci_test_docker.yml`       | Configuration |

Shared helpers (`shared_helpers.py`, `docker_lib.py`, `ci_lib.py`,
`read_yml_config.py`, `setup_logging.py`, `find_py_packages.py`) all live
in `../../shared/` and are imported via a `sys.path` insert at the top of
`ci-oifs-docker.py`.

## Base image

The CI is set up around the official `docker.io/library/gcc:<tag>` images
(Debian-based, GCC pre-installed). To switch to a different base image
family (e.g. plain Ubuntu, or a non-Debian distro), edit `Dockerfile.ci`:
update the `FROM` line and adjust the `apt install` packages — for plain
Ubuntu add `gcc-<N> g++-<N> gfortran-<N>`; for Alpine/RHEL/SUSE swap `apt`
for the matching package manager.

## NORMS convention

The CI relies on the `ifs-test` framework's environment variables (see
`../../ifs-test/README.md`):

- `IFS_TEST_BITIDENTICAL=init` → the framework writes `SAVED_NORMS` per test
- `IFS_TEST_LEGACY=1`          → use the `ifs-grep-norms.pl` legacy flow

## Out of scope

- 3D experiment download/runs (`openifs_expt_url`, `oifs-run`)
- SCM data download/runs (`scm_url`, `callscm`)
- Tolerance-based NORMS comparison (`IFS_TEST_TOLERANCE`) — bit-identical
  `filecmp` only
- CI runner integration (GitHub Actions / GitLab) — only the local scripts
