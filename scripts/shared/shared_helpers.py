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
Generic (non-Docker) helpers shared by the OpenIFS automation drivers
(``bootstrap/docker/create-oifs-docker.py`` and
``ci/docker_ci/ci-oifs-docker.py``) so they interpret ``openifs_source``
values, format timings, slug branch names, clone the OpenIFS repo, copy
local source trees, patch ``oifs-config.edit_me.sh``, and HEAD-check
data URLs the same way.

Docker-specific helpers (image build/pull/validate, Dockerfile rewrite)
live in ``docker_lib.py``.

All functions log via the standard ``logging`` module so the caller controls
where the output goes (file + console set up by ``setup_logging``).
"""

import logging
import os
import re
import shutil
import time
from contextlib import contextmanager
from datetime import timedelta

import git


def format_duration(seconds):
    """Format a duration in seconds as ``H:MM:SS`` for human-readable logs."""
    return str(timedelta(seconds=int(seconds)))


@contextmanager
def timer(description, timings_dict, key):
    """Log a banner around a block of work and record its elapsed time.

    ``timings_dict[key]`` is set to the elapsed seconds when the block exits
    (success or failure), so the driver can print a final timing summary.
    """
    logger = logging.getLogger(__name__)
    logger.info("=" * 70)
    logger.info(description)
    logger.info("=" * 70)
    start_time = time.time()
    try:
        yield
    finally:
        elapsed = time.time() - start_time
        timings_dict[key] = elapsed
        logger.info(f"Completed in {format_duration(elapsed)}")


def shallow_clone(repo_url, clone_dir, branch="main", force=False):
    """Shallow-clone ``branch`` of ``repo_url`` into ``clone_dir``.

    If ``clone_dir`` already exists, the directory is removed when
    ``force=True``; otherwise the clone is skipped (with a warning) so the
    caller can reuse an existing checkout.
    """
    logger = logging.getLogger(__name__)

    if os.path.exists(clone_dir):
        if force:
            logger.info(f"Removing existing directory {clone_dir} (force_reclone=True)")
            shutil.rmtree(clone_dir)
        else:
            logger.warning(f"Directory {clone_dir} already exists")
            logger.info("Skipping clone. Set 'force_reclone: True' in config to override")
            return

    logger.info(f"Cloning {branch} of {repo_url} to {clone_dir}")
    git.Repo.clone_from(repo_url, clone_dir, depth=1, branch=branch)


def check_url_accessible(url, timeout=10):
    """HEAD-request ``url``; return ``(accessible, error_msg)``.

    Used to validate data URLs (ifsdata, climate, rtables, ...) before kicking
    off a long Docker build that would otherwise fail on the first wget.
    """
    import urllib.error
    import urllib.request

    try:
        req = urllib.request.Request(url, method='HEAD')
        with urllib.request.urlopen(req, timeout=timeout) as response:
            if response.status == 200:
                return True, None
            return False, f"HTTP {response.status}"
    except urllib.error.HTTPError as e:
        return False, f"HTTP {e.code}: {e.reason}"
    except urllib.error.URLError as e:
        return False, f"URL Error: {e.reason}"
    except Exception as e:
        return False, f"Error: {str(e)}"


def slug(s):
    """Make a string safe for Docker image tags and filenames.

    Anything outside [A-Za-z0-9._-] is replaced with '-' and leading/trailing
    '-' are trimmed.  Used for branch and directory names that may contain
    '/' or other characters that are awkward in tags or filenames.
    """
    return re.sub(r'[^A-Za-z0-9._-]+', '-', s).strip('-') or 'local'


def resolve_openifs_source(source_value, script_file=None):
    """Resolve an ``openifs_source`` config value to ``(kind, value)``.

    Single source of truth used by all three driver scripts (bootstrap
    docker, bootstrap host, ci docker_ci) so they interpret an
    ``openifs_source`` / ``test_branch`` value identically:

      - empty / None   -> ``('auto', abs_path)``
        The OpenIFS tree that contains ``script_file`` (three directory
        levels up from it).  ``script_file`` must be provided — pass
        ``__file__`` from the calling driver.

      - a directory path -> ``('local', abs_path)``
        An existing local checkout.  ``~`` and ``$VAR`` are expanded.

      - anything else  -> ``('remote', branch_name)``
        A remote branch to clone via ``shallow_clone``.

    Returns ``(kind, value)`` where kind is ``'auto'``, ``'local'``, or
    ``'remote'``.
    """
    source = (source_value or '').strip()
    if not source:
        if script_file is None:
            raise ValueError(
                "openifs_source is empty and script_file was not provided for auto-detection"
            )
        auto_path = os.path.normpath(
            os.path.join(os.path.dirname(os.path.abspath(script_file)), '..', '..', '..')
        )
        return ('auto', auto_path)
    expanded = os.path.expanduser(os.path.expandvars(source))
    if os.path.isdir(expanded):
        return ('local', os.path.normpath(expanded))
    return ('remote', source)


def populate_from_local(local_src, dest_dir, label, force=False):
    """Copy a local OpenIFS checkout into ``dest_dir`` for use as a build source.

    Honours ``force``: when True, an existing ``dest_dir`` is removed and
    re-copied.  When False, an existing ``dest_dir`` is reused as-is (lets
    the user manually populate or rsync between runs).  Skips transient
    artefacts (``.git``, ``build``, ``.cache``, ``.bootstrap``,
    ``__pycache__``, ``openifs-env``).

    If ``local_src`` and ``dest_dir`` resolve to the same path, the copy
    is skipped (the source is already where the build expects it).

    ``label`` is a human-readable tag used only in log messages (e.g.
    ``"control"`` or ``"test"``).
    """
    logger = logging.getLogger(__name__)

    if not os.path.isdir(local_src):
        raise FileNotFoundError(f"Local {label} source {local_src} does not exist")

    if os.path.abspath(local_src) == os.path.abspath(dest_dir):
        logger.info(f"{label.capitalize()} source and build dir are the same ({dest_dir}); skipping copy")
        return

    if os.path.exists(dest_dir) and force:
        logger.info(f"Removing existing {dest_dir} (force_reclone=True)")
        shutil.rmtree(dest_dir)

    if os.path.exists(dest_dir):
        logger.info(f"Using existing {label} source at {dest_dir}")
        return

    logger.info(f"Copying local {label} source {local_src} -> {dest_dir}")
    # '_oifs_docker_ci' is the build/work dir name used by the GitHub Actions
    # CI workflow, which lives inside github.workspace and would otherwise
    # be recursed into as copytree writes into its own subtree.
    shutil.copytree(
        local_src, dest_dir,
        symlinks=True,
        ignore=shutil.ignore_patterns(
            '.git', 'build', '.cache', '.bootstrap',
            '__pycache__', '*.pyc', '*.pyo',
            'openifs-env', '_oifs_docker_ci',
        ),
    )


def move_to_backup(path):
    """Rename ``path`` to ``<path>.backup-<YYYYMMDD-HHMMSS>``.

    Use in place of ``shutil.rmtree`` at user-pickable paths so any
    uncommitted or unpushed work the user happens to have at the path
    survives the operation — fully recoverable with a plain ``mv``.

    Returns the backup path.  If a backup with the same timestamp
    already exists (two runs in the same second, or a stale prior
    backup), a numeric suffix is appended rather than overwriting.
    """
    logger = logging.getLogger(__name__)

    base = f"{path}.backup-{time.strftime('%Y%m%d-%H%M%S')}"
    backup_path = base
    counter = 1
    while os.path.exists(backup_path):
        backup_path = f"{base}.{counter}"
        counter += 1

    os.rename(path, backup_path)
    logger.warning(
        "Moved '%s' aside to '%s'. Any uncommitted or unpushed work is "
        "preserved there.",
        path, backup_path,
    )
    logger.warning(
        "Recover with:  mv '%s' '%s'", backup_path, path,
    )
    logger.warning(
        "Delete old backups once no longer needed:  rm -rf '%s.backup-*'", path,
    )
    return backup_path


def stage_branch_source(config, label, build_dir, script_file):
    """Clone or copy the OpenIFS source for ``label`` into a per-label build dir.

    Single source of truth for the CI driver (docker_ci).  Returns
    ``(clone_dir, tag_suffix)`` where ``clone_dir`` is the absolute path of
    the staged OpenIFS tree and ``tag_suffix`` is a slug safe for image tags
    and filenames (``slug(branch)`` for remote, ``local-<basename>`` for
    local sources).

    For ``label == "test"``, the source is resolved from
    ``config['test_branch']`` via :func:`resolve_openifs_source`:
      - empty / a directory path -> local copy
      - anything else            -> remote clone
    For any other ``label`` (typically ``"control"``), the source is always
    remote — ``config['control_branch']`` cloned from
    ``config['openifs_repo_url']`` — unless ``clone_openifs_control`` is
    False AND a pre-populated tree already exists at the expected path.

    Raises ``FileNotFoundError`` when ``clone_openifs_control`` is False
    and the expected control tree is missing.
    """
    logger = logging.getLogger(__name__)

    branch_build_dir = os.path.join(build_dir, f"build_dir_{label}")
    os.makedirs(branch_build_dir, exist_ok=True)
    clone_dir = os.path.join(branch_build_dir, config['openifs_version'])

    if label == "test":
        kind, value = resolve_openifs_source(config.get('test_branch', ''), script_file)
        if kind != 'remote':
            populate_from_local(value, clone_dir, label,
                                force=config.get('force_reclone', False))
            tag_suffix = f"local-{slug(os.path.basename(os.path.realpath(value)))}"
        else:
            shallow_clone(
                config['openifs_repo_url'], clone_dir,
                branch=value, force=config.get('force_reclone', False),
            )
            tag_suffix = slug(value)
        return clone_dir, tag_suffix

    branch = config['control_branch']
    if config.get('clone_openifs_control', True):
        shallow_clone(
            config['openifs_repo_url'], clone_dir,
            branch=branch, force=config.get('force_reclone', False),
        )
    elif os.path.exists(clone_dir):
        logger.info(f"Using existing {label} source at {clone_dir}")
    else:
        raise FileNotFoundError(
            f"{label.capitalize()} source not found at {clone_dir}. "
            "With clone_openifs_control: False, populate that path yourself "
            "(cp -r / rsync / symlink your working copy), or set "
            "clone_openifs_control: True to clone from the remote."
        )
    tag_suffix = slug(branch)
    return clone_dir, tag_suffix


def patch_oifs_home(oifs_home):
    """Patch ``oifs-config.edit_me.sh`` so ``OIFS_HOME`` points to *oifs_home*.

    Uses a Python regex replace so the edit is portable across macOS and
    Linux (no ``sed -i`` portability issues).  Returns the path to the
    patched config file.
    """
    logger = logging.getLogger(__name__)

    config_file = os.path.join(oifs_home, "oifs-config.edit_me.sh")
    if not os.path.isfile(config_file):
        raise FileNotFoundError(
            f"oifs-config.edit_me.sh not found at {config_file}"
        )

    with open(config_file, 'r') as f:
        content = f.read()

    new_content = re.sub(
        r'^export OIFS_HOME=.*$',
        f'export OIFS_HOME="{oifs_home}"',
        content,
        flags=re.MULTILINE,
    )

    with open(config_file, 'w') as f:
        f.write(new_content)

    logger.info(f"Patched OIFS_HOME in {config_file} -> {oifs_home}")
    return config_file
