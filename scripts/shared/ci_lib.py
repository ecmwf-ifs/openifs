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
CI-specific helpers used by ``ci/docker_ci/ci-oifs-docker.py``.

The driver builds a control branch, builds a test branch, bit-compares
SAVED_NORMS, and writes a self-contained text report. The report-shaped
output, command strings, and summary block live here so they can be
reused if a host-based CI driver is added later.
"""

import os
import subprocess

from shared_helpers import format_duration, resolve_openifs_source, slug


# Standard env prefix that makes openifs-test.sh's framework drop a
# SAVED_NORMS reference file in every test*/ subdirectory.
TEST_ENV_PREFIX = "IFS_TEST_BITIDENTICAL=init IFS_TEST_LEGACY=1"


def _resolve_control_sha(config):
    """Return the 7-char commit SHA of ``control_branch`` on the remote."""
    out = subprocess.check_output(
        ["git", "ls-remote", config['openifs_repo_url'], config['control_branch']],
        text=True,
    ).strip()
    return out.split()[0][:7]


def report_filename(config, script_file):
    """Return the report filename, encoding the two compared branches.

    Format: ``<control>-<sha7>__<test_id>.txt``, where ``<sha7>`` is the
    7-char commit hash of ``control_branch`` resolved via ``git ls-remote``,
    and ``<test_id>`` is either ``slug(test_branch)`` for a remote branch
    or the basename of the resolved local source tree otherwise. Distinct
    comparisons sit alongside each other in ``ci_reports`` without
    overwriting one another.
    """
    control = slug(config['control_branch'])
    sha = _resolve_control_sha(config)
    kind, value = resolve_openifs_source(config.get('test_branch', ''), script_file)
    if kind != 'remote':
        test_id = slug(os.path.basename(os.path.realpath(value)))
    else:
        test_id = slug(value)
    return f"{control}-{sha}__{test_id}.txt"


def control_tarball_name(config, cache_key):
    """Return the filename for the cached control SAVED_NORMS tarball.

    ``cache_key`` distinguishes builds that should not share a cache
    (e.g. ``"gcc14"``): pass whatever string identifies the
    (compiler, version) combination.
    """
    return f"control_saved_norms_{config['openifs_version']}_{cache_key}.tgz"


def build_test_commands(config, source_cmd, output_path):
    """Build the two openifs-test.sh command strings (configure+build, then ctest).

    Split into two stages so the ctest stdout/stderr can be tee'd to
    ``output_path`` in isolation, without dragging in the much louder
    configure/build output.

    ``source_cmd`` is the ``source <oifs-config>`` clause that runs first.
    ``output_path`` is the path that ``tee`` writes to.
    """
    extra_flags = config.get('openifs_test_extra_flags', '').strip()
    cb_cmd = (
        f"{source_cmd} && {TEST_ENV_PREFIX} "
        f"$OIFS_TEST/openifs-test.sh -cb {extra_flags}"
    )
    t_cmd = (
        f"{source_cmd} && {TEST_ENV_PREFIX} "
        f"$OIFS_TEST/openifs-test.sh -t 2>&1 | tee {output_path}"
    )
    return cb_cmd, t_cmd


def write_synthetic_report(report_path, reason):
    """Write a placeholder report when bit-comparison was skipped.

    Used in the control-failure path: there is no bitcompare-generated
    report to copy out, so we synthesise one ourselves so that the rest of
    the report-assembly code (ctest output append, CI summary append) has a
    file to work with.
    """
    os.makedirs(os.path.dirname(report_path), exist_ok=True)
    with open(report_path, "w", encoding="utf-8") as f:
        f.write("=" * 70 + "\n")
        f.write("BIT-COMPARISON SKIPPED\n")
        f.write("=" * 70 + "\n")
        f.write(reason + "\n")
        f.write("\nRESULT: SKIPPED\n")


def append_test_outputs_to_report(report_path, ci_reports, branches):
    """Append captured ctest output for each (label, branch_name) pair to report_path.

    Each section gets a banner so the report is self-documenting. Missing
    files (e.g. control failed before ctest ran) are silently skipped.
    """
    for label, branch_name in branches:
        src = os.path.join(ci_reports, f"openifs_test_output_{label}.txt")
        if not os.path.exists(src):
            continue
        with open(src, encoding="utf-8") as f_in:
            content = f_in.read()
        with open(report_path, "a", encoding="utf-8") as f_out:
            f_out.write("\n")
            f_out.write("=" * 70 + "\n")
            f_out.write(f"CTEST OUTPUT — {label} ({branch_name})\n")
            f_out.write("=" * 70 + "\n")
            f_out.write(content)
            if not content.endswith("\n"):
                f_out.write("\n")


_CONTROL_ANNOTATION = {
    'ok':     'built + tested',
    'reused': 'reused cached NORMS',
    'failed': 'FAILED — bit-comparison skipped',
}


def build_ci_summary(*, control_branch, test_branch, control_status,
                    control_tarball, bit_compare_status, final_status,
                    report_path, timings, total, timing_keys):
    """Build the CI summary lines.

    ``timing_keys`` is the ordered list of timing entries to print. Both
    drivers include ``control-branch``, ``test-branch``, ``norms_compare``;
    docker_ci adds ``base_image`` to the front. Missing keys are skipped.
    """
    annotation = _CONTROL_ANNOTATION[control_status]
    control_norms_line = (
        f"  control NORMS  : {control_tarball} (cached for reuse)"
        if control_status != 'failed'
        else "  control NORMS  : (not produced — control failed)"
    )
    lines = [
        "=" * 70,
        "CI SUMMARY",
        "=" * 70,
        f"  control branch : {control_branch}  [{annotation}]",
        f"  test branch    : {test_branch}",
        f"  bit-comparison : {bit_compare_status}",
        f"  result         : {final_status}",
        f"  report         : {report_path}",
        control_norms_line,
    ]
    for k in timing_keys:
        if k in timings:
            lines.append(f"  {k:<22}: {format_duration(timings[k])}")
    lines.append(f"  total                 : {format_duration(total)}")
    lines.append("=" * 70)
    return lines
