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
    if not out:
        raise ValueError(
            f"Unable to resolve control branch/ref "
            f"'{config['control_branch']}' from remote "
            f"'{config['openifs_repo_url']}': git ls-remote returned no output"
        )
    return out.split()[0][:7]


def resolve_test_branch_label(config, script_file):
    """Return a human-readable label for the test branch in the CI summary.

    For a remote test_branch, returns ``branch (sha7)`` when ``git ls-remote``
    can resolve a SHA, else just the branch name. For a local checkout, queries
    git in the working tree and returns ``branch (sha7)``. Falls back to the
    ``GITHUB_HEAD_REF`` / ``GITHUB_REF_NAME`` env vars when HEAD is detached
    (the usual state after ``actions/checkout``). If git cannot answer at all,
    returns the raw path so the summary still has something to show.
    """
    kind, value = resolve_openifs_source(config.get('test_branch', ''), script_file)
    if kind == 'remote':
        try:
            out = subprocess.check_output(
                ["git", "ls-remote", config['openifs_repo_url'], value],
                text=True, stderr=subprocess.DEVNULL,
            ).strip()
            if out:
                return f"{value} ({out.split()[0][:7]})"
        except subprocess.CalledProcessError:
            pass
        return value

    try:
        sha = subprocess.check_output(
            ["git", "-C", value, "rev-parse", "--short", "HEAD"],
            text=True, stderr=subprocess.DEVNULL,
        ).strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return value

    branch = None
    try:
        head = subprocess.check_output(
            ["git", "-C", value, "rev-parse", "--abbrev-ref", "HEAD"],
            text=True, stderr=subprocess.DEVNULL,
        ).strip()
        if head and head != "HEAD":
            branch = head
    except subprocess.CalledProcessError:
        pass
    if not branch:
        branch = os.environ.get("GITHUB_HEAD_REF") or os.environ.get("GITHUB_REF_NAME")

    return f"{branch} ({sha})" if branch else sha


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


def build_test_commands(config, source_cmd, build_output_path, test_output_path):
    """Build the two openifs-test.sh command strings (configure+build, then ctest).

    Each stage is tee'd to its own file so the captured output survives in
    the uploaded artifact even when one of the stages fails. ``pipefail`` is
    enabled so a non-zero exit from openifs-test.sh propagates through
    ``tee`` — otherwise the shell pipeline would report ``tee``'s exit
    status (0) and silently mask build/test failures.
    """
    extra_flags = config.get('openifs_test_extra_flags', '').strip()
    cb_cmd = (
        f"set -o pipefail; {source_cmd} && {TEST_ENV_PREFIX} "
        f"$OIFS_TEST/openifs-test.sh -cb {extra_flags} 2>&1 | tee {build_output_path}"
    )
    t_cmd = (
        f"set -o pipefail; {source_cmd} && {TEST_ENV_PREFIX} "
        f"$OIFS_TEST/openifs-test.sh -t 2>&1 | tee {test_output_path}"
    )
    return cb_cmd, t_cmd


def write_synthetic_report(report_path, reason):
    """Write a placeholder report when bit-comparison was skipped.

    Used in the control-failure path: there is no bitcompare-generated
    report to copy out, so we synthesise one ourselves so the CI summary can
    still be appended to a report file.
    """
    os.makedirs(os.path.dirname(report_path), exist_ok=True)
    with open(report_path, "w", encoding="utf-8") as f:
        f.write("=" * 70 + "\n")
        f.write("BIT-COMPARISON SKIPPED\n")
        f.write("=" * 70 + "\n")
        f.write(reason + "\n")
        f.write("\nRESULT: SKIPPED\n")


def append_test_outputs_to_report(report_path, ci_reports, branches):
    """Append captured build + ctest output for each (label, branch_name) pair.

    For each label, looks for ``openifs_build_output_<label>.txt`` and
    ``openifs_test_output_<label>.txt`` in ``ci_reports``. Each section
    gets a banner so the report is self-documenting. Missing files (e.g.
    test failed before ctest ran, or build failed before any ctest output
    existed) are silently skipped.
    """
    kinds = [
        ("build",    "BUILD OUTPUT"),
        ("test",     "CTEST OUTPUT"),
        ("lasttest", "CTEST LASTTEST.LOG"),
    ]
    for label, branch_name in branches:
        for kind, banner in kinds:
            src = os.path.join(ci_reports, f"openifs_{kind}_output_{label}.txt")
            if not os.path.exists(src):
                continue
            with open(src, encoding="utf-8") as f_in:
                content = f_in.read()
            with open(report_path, "a", encoding="utf-8") as f_out:
                f_out.write("\n")
                f_out.write("=" * 70 + "\n")
                f_out.write(f"{banner} — {label} ({branch_name})\n")
                f_out.write("=" * 70 + "\n")
                f_out.write(content)
                if not content.endswith("\n"):
                    f_out.write("\n")


_CONTROL_ANNOTATION = {
    'ok':     'built + tested',
    'reused': 'reused cached NORMS',
    'failed': 'FAILED — bit-comparison skipped',
}


def _summary_line(label, value, width=22):
    return f"  {label:<{width}}: {value}"


def build_ci_summary(*, control_branch, test_branch, control_status,
                    control_build_status, control_test_status,
                    test_build_status, test_test_status,
                    control_tarball, bit_compare_status,
                    bit_compare_skip_reason, final_status,
                    report_path, timings, total, timing_keys):
    """Build the CI summary lines.

    ``timing_keys`` is the ordered list of timing entries to print. Both
    drivers include ``control-branch``, ``test-branch``, ``norms_compare``;
    docker_ci adds ``base_image`` to the front. Missing keys are skipped.
    """
    annotation = _CONTROL_ANNOTATION[control_status]
    bit_compare_value = bit_compare_status
    if bit_compare_skip_reason:
        bit_compare_value = f"{bit_compare_status} [{bit_compare_skip_reason}]"

    control_norms = (
        f"{control_tarball} (cached for reuse)"
        if control_status != 'failed'
        else "(not produced — control failed)"
    )
    lines = [
        "=" * 70,
        "CI SUMMARY",
        "=" * 70,
        _summary_line("control branch", f"{control_branch}  [{annotation}]"),
        _summary_line("control Build", control_build_status),
        _summary_line("control Test", control_test_status),
        _summary_line("test branch", test_branch),
        _summary_line("test Build", test_build_status),
        _summary_line("test Test", test_test_status),
        _summary_line("bit-comparison", bit_compare_value),
        _summary_line("bit-comparison result", final_status),
        _summary_line("report", report_path),
        _summary_line("control NORMS", control_norms),
    ]
    for k in timing_keys:
        if k in timings:
            lines.append(_summary_line(k, format_duration(timings[k])))
    lines.append(_summary_line("total", format_duration(total)))
    lines.append("=" * 70)
    return lines
