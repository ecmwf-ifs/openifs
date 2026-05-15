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
Bit-identical comparison of SAVED_NORMS across two OpenIFS branches.

Walks the control tree's ``test*`` subdirs, looks up the same-named subdir
in the test tree, and compares ``SAVED_NORMS`` byte-for-byte with
``filecmp.cmp(..., shallow=False)``. Logs Pass / Fail / Skipped counters.

Invocable standalone via the argparse ``__main__`` block; the CLI exits
non-zero when ``fail_counter > 0``. The ``--report PATH`` flag mirrors all
log output to a file and appends a final ``RESULT: PASS`` / ``RESULT:
FAIL`` line — used by ``ci-oifs-docker.py`` so it can ``docker cp`` a
small text artefact out of the container instead of exporting the full
SAVED_NORMS tree.
"""

import argparse
import filecmp
import logging
import os
import sys


def list_test_dir(base_dir):
    """Return ``{name: full_path}`` for every top-level subdir starting with 'test'."""
    return {
        d: os.path.join(base_dir, d)
        for d in os.listdir(base_dir)
        if os.path.isdir(os.path.join(base_dir, d)) and d.startswith("test")
    }


def test(control_build_test_dir, test_build_test_dir):
    """Compare every test*/SAVED_NORMS pair across the two build trees.

    Walks the control tree's test* subdirs, looks up the same-named subdir in
    the test tree, and compares ``SAVED_NORMS`` byte-for-byte with
    ``filecmp.cmp(..., shallow=False)``. Tests with a missing SAVED_NORMS on
    either side are counted as skipped.

    Returns ``(pass_counter, fail_counter, skipped_counter)``.
    """
    logger = logging.getLogger(__name__)

    logger.info(f"Top level directory for control test output is {control_build_test_dir}")
    logger.info(f"Top level directory for test test output is {test_build_test_dir}")

    logger.info("Create dictionary for control and test runs : key = test name and value = path")
    control_test_dirs = list_test_dir(control_build_test_dir)
    test_test_dirs = list_test_dir(test_build_test_dir)

    pass_counter = 0
    fail_counter = 0
    skipped_counter = 0

    for key, control_path in control_test_dirs.items():

        # If the test isn't present on the test side at all, mirror the
        # upstream "skipped" behaviour rather than treating it as a failure.
        if key not in test_test_dirs:
            logger.warning(f"Test '{key}' missing from test tree - skipping")
            skipped_counter += 1
            continue

        control_saved_norms = os.path.join(control_path, 'SAVED_NORMS')
        test_saved_norms = os.path.join(test_test_dirs[key], 'SAVED_NORMS')

        if os.path.isfile(control_saved_norms) and os.path.isfile(test_saved_norms):

            logger.info(f"Run bit_compare to compare {key} test in control and test")

            bit_compare = filecmp.cmp(control_saved_norms, test_saved_norms, shallow=False)

            if bit_compare:
                logger.info(f"{key} control - test bit comparison : Pass")
                pass_counter += 1
            else:
                logger.warning(f"{key} control - test bit comparison : Fail")
                logger.warning(f"Comparison of {control_path} and {test_test_dirs[key]} SAVED_NORMS, FAILED - check")
                fail_counter += 1
        else:
            logger.warning(f"Either {control_saved_norms} or {test_saved_norms} does not exist")
            logger.warning("Skipping bit comparison test")
            skipped_counter += 1

    logger.warning(f"Bit comparison testing completed, total tests : {len(control_test_dirs)}")
    logger.warning(f"{pass_counter} tests passed, {fail_counter} tests failed and {skipped_counter} tests were skipped ")

    return pass_counter, fail_counter, skipped_counter


def main():

    parser = argparse.ArgumentParser(
        description="Bit-identical comparison of SAVED_NORMS across two OpenIFS branches.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("control_dir",
                        help="Directory containing the control branch's test*/ subdirs (with SAVED_NORMS)")
    parser.add_argument("test_dir",
                        help="Directory containing the test branch's test*/ subdirs (with SAVED_NORMS)")
    parser.add_argument("--report", type=str, default=None,
                        help="Optional path to write a self-contained text report. "
                             "Captures all log output and appends 'RESULT: PASS' or "
                             "'RESULT: FAIL' on the final line.")
    args = parser.parse_args()

    # Console logging is always on. When --report is given, also mirror to a
    # FileHandler so the report file ends up with exactly the same content as
    # stdout. This is what ci-oifs-docker.py copies out of the container.
    handlers = [logging.StreamHandler()]
    if args.report:
        handlers.append(logging.FileHandler(args.report, mode='w', encoding='utf-8'))
    logging.basicConfig(
        level=logging.INFO,
        format='[%(levelname)s] %(name)s.%(funcName)s : %(message)s',
        handlers=handlers,
    )

    _, fail_counter, _ = test(args.control_dir, args.test_dir)

    # Append the machine-readable verdict so callers can grep the report
    # without re-parsing the per-test lines.
    verdict = "PASS" if fail_counter == 0 else "FAIL"
    logging.getLogger(__name__).warning(f"RESULT: {verdict}")

    # Non-zero exit on any failure so this can drive a CI gate directly.
    sys.exit(1 if fail_counter > 0 else 0)


if __name__ == "__main__":
    main()
