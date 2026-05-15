#! /usr/bin/env python3
#
# (C) Copyright 2011- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#
"""
Host-based CI test for OpenIFS — control branch vs test branch.

Mirrors ``ci-oifs-docker.py`` but runs everything directly on the host
(no Docker). Per branch:
  1. Stage the source into ``<openifs_build_host_dir>/build_dir_<label>/<ver>/``
     (clone for remote refs, copy for local sources).
  2. Patch the staged ``oifs-config.edit_me.sh`` so ``OIFS_HOME`` points at
     the staged tree.
  3. Source it and run ``openifs-test.sh -cb`` then ``-t`` with
     ``IFS_TEST_BITIDENTICAL=init IFS_TEST_LEGACY=1`` so the framework
     drops a SAVED_NORMS file in every test*/ subdir. The ctest stage's
    stdout/stderr is tee'd to a separate file in ``ci_reports``.

The control branch is run first and its SAVED_NORMS tree is tarred into
``<control_saved_norms_dir>/control_saved_norms_<openifs_version>_<key>.tgz``.
When ``reuse_control_if_present`` is set and the matching tarball
already exists, the entire control phase is skipped.

Control failure is tolerated, same semantics as the docker driver: the
test phase still runs, the bit-comparison is skipped, a synthetic report
is written, and the script exits 2 (INCONCLUSIVE).

Then the test branch is staged, built, and tested. If the control side
succeeded, the control tarball is extracted to a sibling host directory
and ``openifs_branch_bitcompare.py`` is invoked directly on the host
with the two SAVED_NORMS trees.

Usage:
    python3 ci-oifs-host.py -c config/ci_test_host.yml
"""

import argparse
import logging
import os
import shutil
import subprocess
import sys
import tarfile
import time

# Reuse the shared helpers + CI library.
_SHARED_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "shared")
if _SHARED_DIR not in sys.path:
    sys.path.insert(0, _SHARED_DIR)

import ci_lib
import find_py_packages
import read_yml_config
import setup_logging
import shared_helpers


# The comparator lives in the sibling docker_ci/ directory; it's a pure
# Python script that operates on two directories of NORMS files, so it
# works equally well off-host without modification.
BITCOMPARE_SCRIPT = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)),
                 "..", "docker_ci", "openifs_branch_bitcompare.py")
)

STATUS_PASS = "PASS - Complete and Successful"
STATUS_BUILD_FAILED = "FAILED - Build did not complete"
STATUS_TEST_FAILED = "FAILED - Test did not complete"
STATUS_NO_NORMS = "FAILED - No NORMS produced"
STATUS_NOT_RUN = "SKIPPED - Not run"
STATUS_REUSED_NORMS = "SKIPPED - Reused cached NORMS"


def _fresh_stage_statuses():
    return {"build": STATUS_NOT_RUN, "test": STATUS_NOT_RUN}


def _stages_passed(stage_statuses):
    return (
        stage_statuses["build"] == STATUS_PASS and
        stage_statuses["test"] == STATUS_PASS
    )


def _control_cache_key(config):
    """Cache-key suffix for the control SAVED_NORMS tarball (e.g. ``gcc14``).

    Encodes the compiler family + version so the host-side cache stays
    disjoint from the docker-side cache and from other compiler builds.
    """
    return f"gcc{config['compiler_version']}"


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            "Host CI test for OpenIFS: stage control + test branches, "
            "run openifs-test.sh -cbt in each, and bit-compare SAVED_NORMS "
            "directly on the host."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--config", "-c", type=str, required=True,
                        help="YAML configuration file (see config/ci_test_host.yml)")
    return parser.parse_args()


def _export_lasttest_log(staged_src, ci_reports, label):
    """Copy ctest's per-test detail log from the staged build dir to ``ci_reports``.

    The file is written by ctest at ``<staged_src>/build/Testing/Temporary/
    LastTest.log`` and contains the full stdout/stderr of every test —
    much richer than the wrapper-level ctest output, and the right place
    to look when a test fails.
    """
    logger = logging.getLogger(__name__)
    src = os.path.join(staged_src, "build", "Testing", "Temporary", "LastTest.log")
    dst = os.path.join(ci_reports, f"openifs_lasttest_output_{label}.txt")
    if os.path.exists(src):
        shutil.copyfile(src, dst)
        logger.info(f"Captured {label} LastTest.log -> {dst}")
    else:
        logger.warning(f"No LastTest.log captured for {label} (file missing)")


def run_openifs_tests(staged_src, config, ci_reports, label):
    """Configure + build with ``openifs-test.sh -cb``, then ctest with ``-t``.

    Both stages are tee'd directly to host files in ``ci_reports``, so the
    captured output is preserved even when one stage exits non-zero. The
    ctest LastTest.log is copied in a finally so its detailed per-test
    output also survives a ctest failure.
    """
    logger = logging.getLogger(__name__)
    stage_statuses = _fresh_stage_statuses()

    shared_helpers.patch_oifs_home(staged_src)

    os.makedirs(ci_reports, exist_ok=True)
    build_output = os.path.join(ci_reports, f"openifs_build_output_{label}.txt")
    test_output  = os.path.join(ci_reports, f"openifs_test_output_{label}.txt")

    source_cmd = f"source {staged_src}/oifs-config.edit_me.sh"
    cb_cmd, t_cmd = ci_lib.build_test_commands(config, source_cmd, build_output, test_output)

    try:
        logger.info(f"Configure + build for {label} in {staged_src}")
        try:
            subprocess.run(["bash", "-lc", cb_cmd], cwd=staged_src, check=True)
            stage_statuses["build"] = STATUS_PASS
        except subprocess.CalledProcessError:
            stage_statuses["build"] = STATUS_BUILD_FAILED
            stage_statuses["test"] = "SKIPPED - Build failed"
            return stage_statuses

        logger.info(f"Running ctest for {label} in {staged_src}")
        try:
            subprocess.run(["bash", "-lc", t_cmd], cwd=staged_src, check=True)
            stage_statuses["test"] = STATUS_PASS
        except subprocess.CalledProcessError:
            stage_statuses["test"] = STATUS_TEST_FAILED
    finally:
        _export_lasttest_log(staged_src, ci_reports, label)

    return stage_statuses


def find_saved_norms_root(staged_src):
    """Return the directory whose children are test*/ dirs.

    Walks the staged build tree, finds any ``SAVED_NORMS`` file, and goes
    up one level (the parent of the containing ``test_*`` directory).
    Aborts if nothing is found — that means the bit-identical step never
    wrote references.
    """
    logger = logging.getLogger(__name__)
    for root, _dirs, files in os.walk(staged_src):
        if "SAVED_NORMS" in files:
            test_root = os.path.dirname(root)
            logger.info(f"SAVED_NORMS root at {test_root}")
            return test_root
    raise FileNotFoundError(
        f"No SAVED_NORMS found under {staged_src} - tests did not produce reference NORMS"
    )


def export_control_norms(test_root, control_tarball):
    """Tar the control SAVED_NORMS tree at ``test_root`` into ``control_tarball``.

    Counterpart to the docker driver's ``docker exec tar … && docker cp``
    sequence: here both steps are local, so we just walk and pack directly.
    """
    logger = logging.getLogger(__name__)
    os.makedirs(os.path.dirname(control_tarball), exist_ok=True)
    logger.info(f"Bundling control SAVED_NORMS -> {control_tarball}")
    with tarfile.open(control_tarball, "w:gz") as tar:
        for name in sorted(os.listdir(test_root)):
            tar.add(os.path.join(test_root, name), arcname=name)


def run_control_phase(config, build_dir, control_tarball, ci_reports):
    """Stage the control branch, build, run tests, export SAVED_NORMS tarball.

    Tolerates failure: if any step raises (clone error or ctest failure),
    we log it and return ``'failed'`` — the caller then skips the
    bit-comparison and reports INCONCLUSIVE.
    """
    logger = logging.getLogger(__name__)

    if config.get('reuse_control_if_present', False) and os.path.exists(control_tarball):
        logger.info("=" * 70)
        logger.info(f"Reusing existing control tarball: {control_tarball}")
        logger.info("(reuse_control_if_present=True; delete the tarball or set the flag")
        logger.info(" to False to force a fresh control run)")
        logger.info("=" * 70)
        return 'reused', {"build": STATUS_REUSED_NORMS, "test": STATUS_REUSED_NORMS}

    stage_statuses = _fresh_stage_statuses()
    try:
        clone_dir, _ = shared_helpers.stage_branch_source(
            config, "control", build_dir, __file__,
        )
        stage_statuses = run_openifs_tests(clone_dir, config, ci_reports, "control")
        if not _stages_passed(stage_statuses):
            return 'failed', stage_statuses
        test_root = find_saved_norms_root(clone_dir)
        export_control_norms(test_root, control_tarball)
        return 'ok', stage_statuses
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        logger.error(f"Control phase FAILED: {e}")
        if stage_statuses["build"] == STATUS_NOT_RUN:
            stage_statuses["build"] = STATUS_BUILD_FAILED
        if _stages_passed(stage_statuses):
            stage_statuses["test"] = STATUS_NO_NORMS
        return 'failed', stage_statuses


def run_test_phase(config, build_dir, ci_reports):
    """Stage the test branch, build, run tests, locate SAVED_NORMS root.

    Symmetric to ``run_control_phase`` but does NOT tolerate failures — a
    test-side build/ctest failure is a hard error and aborts the script.
    """
    clone_dir, _ = shared_helpers.stage_branch_source(
        config, "test", build_dir, __file__,
    )
    stage_statuses = run_openifs_tests(clone_dir, config, ci_reports, "test")
    if not _stages_passed(stage_statuses):
        return None, stage_statuses
    try:
        test_root = find_saved_norms_root(clone_dir)
    except FileNotFoundError:
        stage_statuses["test"] = STATUS_NO_NORMS
        return None, stage_statuses
    return test_root, stage_statuses


def compare_norms(test_root, control_tarball, report_path, build_dir):
    """Extract the control tarball and run the comparator on host paths.

    Returns True iff the comparator exited 0.
    """
    logger = logging.getLogger(__name__)

    control_extract_dir = os.path.join(build_dir, "control_saved_norms_extracted")
    if os.path.exists(control_extract_dir):
        shutil.rmtree(control_extract_dir)
    os.makedirs(control_extract_dir, exist_ok=True)
    logger.info(f"Extracting {control_tarball} -> {control_extract_dir}")
    with tarfile.open(control_tarball, "r:gz") as tar:
        tar.extractall(control_extract_dir)

    os.makedirs(os.path.dirname(report_path), exist_ok=True)
    logger.info(f"Running {BITCOMPARE_SCRIPT}")
    # Non-zero exit means "tests disagree" — that's a normal CI outcome,
    # so we do NOT pass check=True.
    result = subprocess.run(
        ["python3", BITCOMPARE_SCRIPT,
         control_extract_dir, test_root,
         "--report", report_path],
    )
    return result.returncode == 0


def main():
    script_start = time.time()
    timings = {}

    cli_args = parse_arguments()

    find_py_packages.main(["yaml"])

    config = read_yml_config.main(cli_args.config)

    build_dir = config['openifs_build_host_dir']
    os.makedirs(build_dir, exist_ok=True)

    log_dir = os.path.join(build_dir, "host_ci_logfiles")
    os.makedirs(log_dir, exist_ok=True)
    log_path = os.path.join(
        log_dir,
        f"log_ci_{config['openifs_version']}_gcc{config['compiler_version']}.log",
    )
    setup_logging.main(log_path)
    logger = logging.getLogger(__name__)

    ci_reports = config['ci_reports']
    control_dir = config['control_saved_norms_dir']
    control_tarball = os.path.join(
        control_dir,
        ci_lib.control_tarball_name(config, _control_cache_key(config)),
    )
    report_path = os.path.join(ci_reports, ci_lib.report_filename(config, __file__))

    # --- Control phase (skipped if tarball exists and reuse=True) ----------
    with shared_helpers.timer(f"Control phase ({config['control_branch']})", timings, 'control-branch'):
        control_status, control_stage_statuses = run_control_phase(
            config, build_dir, control_tarball, ci_reports,
        )

    # --- Test phase: always runs, even if control failed -------------------
    # A build/ctest failure here is a real CI failure, but we catch it so
    # the captured outputs still land in the uploaded artifact.
    test_branch = config['test_branch']
    test_label = test_branch or "auto-resolved local source"
    test_root = None
    test_stage_statuses = _fresh_stage_statuses()
    test_phase_failed = False
    try:
        with shared_helpers.timer(f"Test phase ({test_label})", timings, 'test-branch'):
            test_root, test_stage_statuses = run_test_phase(config, build_dir, ci_reports)
            test_phase_failed = not _stages_passed(test_stage_statuses)
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        logger.error(f"Test phase FAILED: {e}")
        test_phase_failed = True
        test_stage_statuses["build"] = STATUS_BUILD_FAILED

    # --- Bit-comparison: only if both phases produced SAVED_NORMS ----------
    if test_phase_failed:
        ci_lib.write_synthetic_report(
            report_path,
            "Test phase FAILED during configure+build or ctest. "
            "See the uploaded BUILD OUTPUT and CTEST OUTPUT artifacts for the cause.",
        )
        bit_compare_status = 'SKIPPED'
        bit_compare_skip_reason = 'No test NORMS'
        timings['norms_compare'] = 0
    elif control_status == 'failed':
        logger.warning("Skipping bit-comparison — control phase did not produce SAVED_NORMS")
        ci_lib.write_synthetic_report(
            report_path,
            "Control phase FAILED — no SAVED_NORMS to compare against. "
            "Test phase ran to completion; see uploaded artifacts for details.",
        )
        bit_compare_status = 'SKIPPED'
        bit_compare_skip_reason = 'No control NORMS'
        timings['norms_compare'] = 0
    else:
        with shared_helpers.timer("NORMS comparison (host)", timings, 'norms_compare'):
            passed = compare_norms(test_root, control_tarball, report_path, build_dir)
        bit_compare_status = 'PASS' if passed else 'FAIL'
        bit_compare_skip_reason = None

    total = time.time() - script_start

    # Final result classification:
    #   PASS         control produced NORMS, test bit-matched control
    #   FAIL         control produced NORMS, test bit-DIFFERED, OR test phase failed
    #   INCONCLUSIVE control failed; comparison skipped
    if test_phase_failed:
        final_status, exit_code = 'FAIL', 1
    elif control_status == 'failed':
        final_status, exit_code = 'INCONCLUSIVE', 2
    elif bit_compare_status == 'PASS':
        final_status, exit_code = 'PASS', 0
    else:
        final_status, exit_code = 'FAIL', 1

    summary_lines = ci_lib.build_ci_summary(
        control_branch=config['control_branch'],
        test_branch=ci_lib.resolve_test_branch_label(config, __file__),
        control_status=control_status,
        control_build_status=control_stage_statuses["build"],
        control_test_status=control_stage_statuses["test"],
        test_build_status=test_stage_statuses["build"],
        test_test_status=test_stage_statuses["test"],
        control_tarball=control_tarball,
        bit_compare_status=bit_compare_status,
        bit_compare_skip_reason=bit_compare_skip_reason,
        final_status=final_status,
        report_path=report_path,
        timings=timings,
        total=total,
        timing_keys=('control-branch', 'test-branch', 'norms_compare'),
    )

    for line in summary_lines:
        logger.info(line)

    try:
        with open(report_path, "a", encoding="utf-8") as f:
            f.write("\n")
            for line in summary_lines:
                f.write(line + "\n")
    except OSError as e:
        logger.warning(f"Could not append CI summary to {report_path}: {e}")

    sys.exit(exit_code)


if __name__ == "__main__":
    main()
