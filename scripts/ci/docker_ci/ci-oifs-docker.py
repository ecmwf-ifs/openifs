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
Docker-based CI test for OpenIFS — control branch vs test branch.

Per branch:
  1. Source the branch into a per-label build directory (clone or local copy)
  2. Build a Docker image from Dockerfile.ci (no SCM / experiment data)
  3. Configure + build with ``openifs-test.sh -cb``, then run ctest with
     ``-t``. Both calls set IFS_TEST_BITIDENTICAL=init IFS_TEST_LEGACY=1 so
     the framework drops a SAVED_NORMS file in every test*/ subdir. The
    ctest stage's stdout/stderr is tee'd to a file inside the container
    and ``docker cp``'d out as a separate uploaded artifact.

The control branch is run first and its SAVED_NORMS tree is bundled out as
``<control_saved_norms_dir>/control_saved_norms_<openifs_version>_gcc<base_image>.tgz``
(intermediate blob — kept for debugging, never read by host code). The
name is keyed on the OpenIFS version and GCC base image so different
(version, compiler) combinations cache independently. When
``reuse_control_if_present`` is set in the config and the matching tarball
already exists, the entire control phase is skipped.

Control failure is tolerated. If the control build/ctest fails (e.g. main
doesn't compile and the test branch is the fix), the script:
  - logs the failure, cleans up the control container
  - continues to the test phase (which is still strict — test failures abort)
  - skips the bit-comparison
  - writes a synthetic report explaining the situation
  - exits with code 2 (INCONCLUSIVE) so a CI gate can distinguish from a
    clean PASS (0) or clean FAIL (1)

Then the test branch is built and run. If the control side succeeded:
  4. ``docker cp`` the control tarball INTO the test container
  5. ``docker cp`` openifs_branch_bitcompare.py INTO the test container
  6. ``docker exec`` the comparator inside the test container with --report —
     it untars the control tree, compares each test*/SAVED_NORMS pair, and
     writes a self-contained text report
  7. ``docker cp`` ONLY the report out to
     ``<ci_reports>/<control>-<sha7>__<test|dir>.txt``

The host never reads NORMS data. A final CI summary is appended to the report
file. Captured build and ctest outputs remain separate files in the uploaded
artifact. The script's exit code drives PASS / FAIL / INCONCLUSIVE.

Usage:
    python3 ci-oifs-docker.py -c config/ci_test_docker.yml
"""

import argparse
import logging
import os
import shutil
import subprocess
import sys
import time

# Generic helpers (YAML loader, logger setup, module checker, shared_helpers)
# and Docker-specific helpers (docker_lib) all live in scripts/shared/.
_SHARED_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "shared")
if _SHARED_DIR not in sys.path:
    sys.path.insert(0, _SHARED_DIR)

import ci_lib
import docker_lib
import find_py_packages
import read_yml_config
import setup_logging
import shared_helpers


# Filename of this script's local comparator — copied INTO the test container.
BITCOMPARE_SCRIPT = "openifs_branch_bitcompare.py"

# Inside-container paths used by the containerised comparison step.
INCONTAINER_CONTROL_TGZ = "/tmp/control_saved_norms.tgz"
INCONTAINER_CONTROL_DIR = "/tmp/control_saved_norms"
INCONTAINER_BITCOMPARE = f"/tmp/{BITCOMPARE_SCRIPT}"
INCONTAINER_REPORT = "/tmp/norms_report.txt"

# Inside-container paths that capture each stage's stdout/stderr. Copied
# out per-branch by export_build_output() / export_test_output() so the
# uploaded artifact includes them even when the stage failed.
INCONTAINER_BUILD_OUTPUT = "/tmp/openifs_build_output.txt"
INCONTAINER_TEST_OUTPUT  = "/tmp/openifs_test_output.txt"

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

    Keys the cache on the GCC base image so that different compilers don't
    silently share the previous run's NORMS.
    """
    return f"gcc{config['base_docker_image']}"


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            "CI test for OpenIFS: build control + test branches in Docker, "
            "run openifs-test.sh -cbt in each, and bit-compare SAVED_NORMS "
            "INSIDE the test container. Only a small text report is exported "
            "to the host."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--config", "-c", type=str, required=True,
                        help="YAML configuration file (see config/ci_test_docker.yml)")
    return parser.parse_args()


def ensure_base_image(config):
    """Validate + pull the GCC base image. Mirrors create-oifs-docker.py."""
    logger = logging.getLogger(__name__)
    base_image = f"gcc:{config['base_docker_image']}"

    logger.info(f"Validating base Docker image {base_image}...")
    if not docker_lib.is_official_docker_image(base_image):
        logger.error(f"Security check failed: '{base_image}' is not an approved official image")
        sys.exit(1)

    if not docker_lib.check_docker_image_exists(base_image):
        logger.info(f"Base image {base_image} not present locally, pulling...")
        if not docker_lib.pull_docker_image(base_image):
            logger.error(f"Failed to pull {base_image}")
            sys.exit(1)


def build_branch_image(config, label, branch, build_dir):
    """Stage the source for ``label`` into a per-label build dir and build the image.

    Source staging (clone vs. local copy, force_reclone, tag-suffix
    derivation) is delegated to ``shared_helpers.stage_branch_source`` —
    same routine used by ``ci-oifs-host.py``. Each branch needs its own
    Docker build dir because Dockerfile.ci uses
    ``COPY ${OPENIFS_DIR} ...`` and ``docker build`` ships the entire
    build dir to the daemon.

    Returns the image tag.
    """
    logger = logging.getLogger(__name__)

    try:
        clone_dir, tag_suffix = shared_helpers.stage_branch_source(
            config, label, build_dir, __file__,
        )
    except FileNotFoundError as e:
        logger.error(str(e))
        sys.exit(1)

    branch_build_dir = os.path.dirname(clone_dir)
    branch_dockerfile = os.path.join(
        branch_build_dir,
        f"Dockerfile_ci_{config['openifs_version']}_{config['base_docker_image']}",
    )

    shutil.copyfile(config['docker_template'], branch_dockerfile)
    branch_config = dict(
        config,
        openifs_branch=tag_suffix,
        include_openifs_data_downloads=False,
    )
    docker_lib.modify_dockerfile(branch_dockerfile, branch_config)

    image_tag = (
        f"openifs-{config['openifs_version']}-gcc{config['base_docker_image']}"
        f":ci-{label}-{tag_suffix}"
    )
    logger.info(f"Building image {image_tag} for {label} branch '{branch}'...")
    docker_lib.build_docker_image(
        branch_dockerfile,
        image_tag,
        branch_build_dir,
        no_cache=config.get('force_rebuild', False),
    )
    return image_tag


def start_container(label, image_tag, config):
    """Start a detached bash container; return its name.

    By default, removes any pre-existing container with the same name first
    (e.g. from a prior failed run). This behaviour is controlled by
    ``remove_existing_container_before_run`` in the YAML and defaults to True.
    """
    logger = logging.getLogger(__name__)

    container = f"oifs-ci-{label}"
    if config.get('remove_existing_container_before_run', True):
        subprocess.run(
            ["docker", "rm", "-f", container],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        logger.info(f"Removed any pre-existing container named '{container}' before run")

    logger.info(f"Starting container '{container}' from {image_tag}")
    subprocess.run(
        ["docker", "run", "-dit", "--name", container, image_tag, "/bin/bash"],
        check=True,
    )
    return container


def run_openifs_tests(container, config, ci_reports, label):
    """Configure + build with ``openifs-test.sh -cb``, then run ctest with ``-t``.

    Both stages are tee'd to per-stage files inside the container. The
    finally block always copies those files out to ``ci_reports`` — even
    when one of the stages exited non-zero — so the failure cause is
    visible in the uploaded artifact.
    """
    logger = logging.getLogger(__name__)
    stage_statuses = _fresh_stage_statuses()

    src = f"source ~/{config['openifs_version']}/oifs-config.edit_me.sh"
    cb_cmd, t_cmd = ci_lib.build_test_commands(
        config, src, INCONTAINER_BUILD_OUTPUT, INCONTAINER_TEST_OUTPUT,
    )

    try:
        logger.info(f"Configure + build in '{container}'")
        try:
            subprocess.run(
                ["docker", "exec", container, "bash", "-lc", cb_cmd],
                check=True,
            )
            stage_statuses["build"] = STATUS_PASS
        except subprocess.CalledProcessError:
            stage_statuses["build"] = STATUS_BUILD_FAILED
            stage_statuses["test"] = "SKIPPED - Build failed"
            return stage_statuses

        logger.info(f"Running ctest in '{container}'")
        try:
            subprocess.run(
                ["docker", "exec", container, "bash", "-lc", t_cmd],
                check=True,
            )
            stage_statuses["test"] = STATUS_PASS
        except subprocess.CalledProcessError:
            stage_statuses["test"] = STATUS_TEST_FAILED
    finally:
        export_build_output(container, ci_reports, label)
        export_test_output(container, ci_reports, label)
        export_lasttest_log(container, config, ci_reports, label)

    return stage_statuses


def _export_stage_output(container, ci_reports, label, in_container_path, kind):
    """Copy a per-stage output file from ``container`` into ``ci_reports``.

    ``kind`` is ``"build"`` or ``"test"``. The host-side filename is
    ``openifs_<kind>_output_<label>.txt``. Returns the host path on
    success, or None if the in-container file is missing (e.g. the stage
    never ran because an earlier stage failed).
    """
    logger = logging.getLogger(__name__)

    os.makedirs(ci_reports, exist_ok=True)
    host_path = os.path.join(ci_reports, f"openifs_{kind}_output_{label}.txt")
    try:
        subprocess.run(
            ["docker", "cp", f"{container}:{in_container_path}", host_path],
            check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )
        logger.info(f"Captured {label} {kind} output -> {host_path}")
        return host_path
    except subprocess.CalledProcessError:
        logger.warning(f"No {kind} output captured for {label} (file missing in container)")
        return None


def export_build_output(container, ci_reports, label):
    """Copy the captured configure+build output out of ``container``."""
    return _export_stage_output(container, ci_reports, label,
                                INCONTAINER_BUILD_OUTPUT, "build")


def export_test_output(container, ci_reports, label):
    """Copy the captured ctest output out of ``container``."""
    return _export_stage_output(container, ci_reports, label,
                                INCONTAINER_TEST_OUTPUT, "test")


def export_lasttest_log(container, config, ci_reports, label):
    """Copy ctest's per-test detail log out of ``container``.

    Located at ``OIFS_HOME/build/Testing/Temporary/LastTest.log`` inside the
    container, where ``OIFS_HOME = /home/openifs/<openifs_version>`` per
    Dockerfile.ci. This file is ctest-generated and contains the full
    stdout/stderr of every individual test — much richer than the
    wrapper-level ctest output, and the right place to look when a test
    fails.
    """
    in_container_path = (
        f"/home/openifs/{config['openifs_version']}"
        f"/build/Testing/Temporary/LastTest.log"
    )
    return _export_stage_output(container, ci_reports, label,
                                in_container_path, "lasttest")


def find_saved_norms_root(container):
    """Return the in-container directory whose children are test*/ dirs.

    Located by finding any SAVED_NORMS file and going up two levels
    (test*/SAVED_NORMS -> test_dir -> shared parent). Aborts the run if
    nothing is found, since that means the bit-identical step never wrote
    references.
    """
    logger = logging.getLogger(__name__)

    locate_cmd = (
        "set -e && cd $OIFS_HOME && "
        "TEST_ROOT=$(find . -type f -name SAVED_NORMS -print -quit | xargs -n1 dirname | xargs -n1 dirname) && "
        "[ -n \"$TEST_ROOT\" ] || { echo 'No SAVED_NORMS found - tests did not produce reference NORMS' >&2; exit 1; } && "
        "(cd \"$TEST_ROOT\" && pwd)"
    )
    result = subprocess.run(
        ["docker", "exec", container, "bash", "-lc", locate_cmd],
        check=True, capture_output=True, text=True,
    )
    test_root = result.stdout.strip().splitlines()[-1]
    logger.info(f"SAVED_NORMS root inside '{container}': {test_root}")
    return test_root


def export_control_norms(container, test_root, control_tarball):
    """Tar the control SAVED_NORMS tree and copy the tarball out to the host.

    This is the only intermediate that touches the host filesystem — it's an
    opaque blob that gets shipped INTO the test container later for
    comparison. The host never reads it. Its presence on disk is also the
    "successful control run" signal for ``reuse_control_if_present``.
    """
    logger = logging.getLogger(__name__)

    bundle_cmd = f"tar czf /tmp/saved_norms.tgz -C \"{test_root}\" ."
    subprocess.run(
        ["docker", "exec", container, "bash", "-lc", bundle_cmd],
        check=True,
    )

    os.makedirs(os.path.dirname(control_tarball), exist_ok=True)
    logger.info(f"Copying control SAVED_NORMS tarball to {control_tarball}")
    subprocess.run(
        ["docker", "cp", f"{container}:/tmp/saved_norms.tgz", control_tarball],
        check=True,
    )


def remove_container_if_requested(container, config):
    """Honour ``remove_test_container`` from the YAML."""
    logger = logging.getLogger(__name__)
    if config.get('remove_test_container', True):
        subprocess.run(["docker", "rm", "-f", container], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        logger.info(f"Removed container '{container}'")
    else:
        logger.info(f"Container '{container}' left running for inspection")


def run_control_phase(config, build_dir, control_tarball, ci_reports):
    """Build the control branch, run tests, export SAVED_NORMS tarball.

    Tolerates failure: if the clone/build/ctest steps raise (a docker exec
    returned non-zero), we log it, clean up the container, and return
    ``'failed'`` — the caller should then skip the bit-comparison and run
    in INCONCLUSIVE mode rather than aborting.

    Returns one of:
            - ``('reused', statuses)`` — short-circuited via cache
            - ``('ok', statuses)``     — built + tested + exported successfully
            - ``('failed', statuses)`` — a step failed; tarball was NOT written
    """
    logger = logging.getLogger(__name__)

    if config.get('reuse_control_if_present', False) and os.path.exists(control_tarball):
        logger.info("=" * 70)
        logger.info(f"Reusing existing control tarball: {control_tarball}")
        logger.info("(reuse_control_if_present=True; delete the tarball or set the flag")
        logger.info(" to False to force a fresh control run)")
        logger.info("=" * 70)
        return 'reused', {"build": STATUS_REUSED_NORMS, "test": STATUS_REUSED_NORMS}

    branch = config['control_branch']
    container = None
    stage_statuses = _fresh_stage_statuses()
    try:
        image_tag = build_branch_image(config, "control", branch, build_dir)
        container = start_container("control", image_tag, config)
        stage_statuses = run_openifs_tests(container, config, ci_reports, "control")
        if not _stages_passed(stage_statuses):
            return 'failed', stage_statuses
        test_root = find_saved_norms_root(container)
        export_control_norms(container, test_root, control_tarball)
        return 'ok', stage_statuses
    except subprocess.CalledProcessError as e:
        logger.error(f"Control phase FAILED: {e}")
        # run_openifs_tests' finally block already exported whatever build /
        # ctest output exists, even on failure. Nothing more to do here.
        if stage_statuses["build"] == STATUS_NOT_RUN:
            stage_statuses["build"] = STATUS_BUILD_FAILED
        if _stages_passed(stage_statuses):
            stage_statuses["test"] = STATUS_NO_NORMS
        return 'failed', stage_statuses
    finally:
        if container is not None:
            remove_container_if_requested(container, config)


def run_test_phase(config, build_dir, ci_reports):
    """Build the test branch, run tests, capture ctest output, locate SAVED_NORMS.

    Symmetric to ``run_control_phase`` but does NOT tolerate failures — a
    test-side build/ctest failure is a hard error and aborts the script.
    The container is deliberately left RUNNING so the subsequent
    in-container bit-comparison step can copy the control tarball and the
    comparator script into it. The caller is responsible for removing the
    container afterwards.

    Returns ``(container, test_root, statuses)``.
    """
    branch = config['test_branch']
    image_tag = build_branch_image(config, "test", branch, build_dir)
    container = start_container("test", image_tag, config)
    stage_statuses = run_openifs_tests(container, config, ci_reports, "test")
    if not _stages_passed(stage_statuses):
        return container, None, stage_statuses
    try:
        test_root = find_saved_norms_root(container)
    except subprocess.CalledProcessError:
        stage_statuses["test"] = STATUS_NO_NORMS
        return container, None, stage_statuses
    return container, test_root, stage_statuses


def compare_norms_in_container(test_container, test_root, control_tarball, report_path):
    """Run the bit-comparison INSIDE ``test_container`` and copy the report out.

    Ships the control tarball + the comparator script into the test
    container's /tmp, runs ``openifs_branch_bitcompare.py --report``, and
    copies the small text report back to the host at the caller-provided
    ``report_path``.

    Returns True iff the comparator exited 0.
    """
    logger = logging.getLogger(__name__)

    # Ship the control SAVED_NORMS blob + the comparator into the test
    # container. Both go to /tmp because that's writable for the openifs user.
    logger.info(f"Copying {control_tarball} into '{test_container}':{INCONTAINER_CONTROL_TGZ}")
    subprocess.run(
        ["docker", "cp", control_tarball, f"{test_container}:{INCONTAINER_CONTROL_TGZ}"],
        check=True,
    )

    bitcompare_src = os.path.join(os.path.dirname(os.path.abspath(__file__)), BITCOMPARE_SCRIPT)
    logger.info(f"Copying {bitcompare_src} into '{test_container}':{INCONTAINER_BITCOMPARE}")
    subprocess.run(
        ["docker", "cp", bitcompare_src, f"{test_container}:{INCONTAINER_BITCOMPARE}"],
        check=True,
    )

    # Run the comparator inside the test container. We DON'T pass check=True
    # — a non-zero exit code is the comparator's way of saying "tests
    # disagree", which is a normal CI outcome, not an orchestration error.
    compare_cmd = (
        f"set -e && "
        f"rm -rf {INCONTAINER_CONTROL_DIR} && mkdir -p {INCONTAINER_CONTROL_DIR} && "
        f"tar xzf {INCONTAINER_CONTROL_TGZ} -C {INCONTAINER_CONTROL_DIR} && "
        # The bitcompare script returns non-zero on FAIL — let that bubble up.
        f"set +e && python3 {INCONTAINER_BITCOMPARE} "
        f"{INCONTAINER_CONTROL_DIR} \"{test_root}\" --report {INCONTAINER_REPORT}"
    )
    logger.info(f"Running comparator in '{test_container}'")
    result = subprocess.run(
        ["docker", "exec", test_container, "bash", "-lc", compare_cmd],
    )
    passed = (result.returncode == 0)

    # Always pull the report — even on failure, the user wants to see why.
    os.makedirs(os.path.dirname(report_path), exist_ok=True)
    logger.info(f"Copying report to {report_path}")
    subprocess.run(
        ["docker", "cp", f"{test_container}:{INCONTAINER_REPORT}", report_path],
        check=True,
    )

    return passed


def main():
    script_start = time.time()
    timings = {}

    cli_args = parse_arguments()

    # pyyaml is lazy-imported by read_yml_config, so check it here for an
    # early-fail with a friendly message. (gitpython is already enforced by
    # the top-level `import shared_helpers`, which imports git itself.)
    find_py_packages.main(["yaml"])

    config = read_yml_config.main(cli_args.config)

    build_dir = config['openifs_build_docker_dir']
    os.makedirs(build_dir, exist_ok=True)

    log_dir = os.path.join(build_dir, "docker_ci_logfiles")
    os.makedirs(log_dir, exist_ok=True)
    log_path = os.path.join(
        log_dir,
        f"log_ci_{config['openifs_version']}_{config['base_docker_image']}.log",
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

    # Resolve docker_template relative to this file so the YAML can use a
    # short "./Dockerfile.ci" without depending on the user's CWD.
    if not os.path.isabs(config['docker_template']):
        config['docker_template'] = os.path.normpath(
            os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         config['docker_template'])
        )

    # Pre-flight: base image security check + pull
    with shared_helpers.timer("Base image validation", timings, 'base_image'):
        ensure_base_image(config)

    # --- Control phase (skipped if tarball already exists and reuse=True) ------
    # Returns 'ok', 'reused', or 'failed'. A 'failed' control still lets the
    # rest of the script run — we just skip the bit-comparison and report
    # INCONCLUSIVE. Useful when the test branch is the fix for a broken main.
    with shared_helpers.timer(f"Control phase ({config['control_branch']})", timings, 'control-branch'):
        control_status, control_stage_statuses = run_control_phase(
            config, build_dir, control_tarball, ci_reports,
        )

    # --- Test phase: always runs, even if the control phase failed ------------
    # A build/ctest failure here is a real CI failure, but we catch it so
    # the captured outputs still land in the uploaded artifact.
    test_branch = config['test_branch']
    test_label = test_branch or "auto-resolved local source"
    test_container = None
    test_root = None
    test_stage_statuses = _fresh_stage_statuses()
    test_phase_failed = False
    try:
        with shared_helpers.timer(f"Test phase ({test_label})", timings, 'test-branch'):
            test_container, test_root, test_stage_statuses = run_test_phase(
                config, build_dir, ci_reports,
            )
            test_phase_failed = not _stages_passed(test_stage_statuses)
    except subprocess.CalledProcessError as e:
        logger.error(f"Test phase FAILED: {e}")
        test_phase_failed = True
        test_stage_statuses["build"] = STATUS_BUILD_FAILED

    # --- Bit-comparison: only if both phases produced SAVED_NORMS -------------
    if test_phase_failed:
        ci_lib.write_synthetic_report(
            report_path,
            "Test phase FAILED during configure+build or ctest. "
            "See the uploaded BUILD OUTPUT and CTEST OUTPUT artifacts for the cause.",
        )
        bit_compare_status = 'SKIPPED'
        bit_compare_skip_reason = 'No test NORMS'
        timings['norms_compare'] = 0
        if test_container is not None:
            remove_container_if_requested(test_container, config)
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
        remove_container_if_requested(test_container, config)
    else:
        with shared_helpers.timer("NORMS comparison (in container)", timings, 'norms_compare'):
            passed = compare_norms_in_container(
                test_container, test_root, control_tarball, report_path,
            )
            remove_container_if_requested(test_container, config)
        bit_compare_status = 'PASS' if passed else 'FAIL'
        bit_compare_skip_reason = None

    total = time.time() - script_start

    # Final result classification:
    #   PASS         control produced NORMS, test bit-matched control
    #   FAIL         control produced NORMS, test bit-DIFFERED, OR test phase failed
    #   INCONCLUSIVE control failed; comparison skipped (test ran fine on its own)
    if test_phase_failed:
        final_status, exit_code = 'FAIL', 1
    elif control_status == 'failed':
        final_status, exit_code = 'INCONCLUSIVE', 2
    elif bit_compare_status == 'PASS':
        final_status, exit_code = 'PASS', 0
    else:
        final_status, exit_code = 'FAIL', 1

    # Build the summary once so we can both log it and append it to the report.
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
        timing_keys=('base_image', 'control-branch', 'test-branch', 'norms_compare'),
    )

    for line in summary_lines:
        logger.info(line)

    # Append the summary to the report so it is fully self-contained. Earlier
    # content (bitcompare body or synthetic header) is preserved.
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
