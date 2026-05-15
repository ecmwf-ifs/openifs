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
Docker-specific helpers used by the OpenIFS Docker drivers.

Used by ``scripts/bootstrap/docker/create-oifs-docker.py`` (single-branch
builder) and ``scripts/ci/docker_ci/ci-oifs-docker.py`` (control-vs-test CI):

  - Docker base-image validation / pull / existence check
  - Dockerfile template rewrite (substitute ARGs and base image)
  - Docker image build

Generic helpers (timing, slugging, git clone, URL HEAD-check, OpenIFS
source resolution) live in ``shared_helpers.py``.

All functions log via the standard ``logging`` module so the caller controls
where the output goes (file + console set up by ``setup_logging``).
"""

import logging
import subprocess
import sys

from shared_helpers import check_url_accessible


# Whitelist used by is_official_docker_image() — only base images from these
# official Docker library names are accepted as the build base image.
ALLOWED_OFFICIAL_IMAGES = [
    'gcc',
    'ubuntu',
    'debian',
]


def is_official_docker_image(image_name):
    """Return True if ``image_name`` is from the Docker Hub official library.

    Accepts the bare form ``gcc:14``, the explicit ``library/gcc:14`` form,
    and the fully-qualified ``docker.io/library/gcc:14`` form. Anything else
    (user/org images, third-party registries) is rejected — the driver uses
    this as a security gate before pulling the base image.
    """
    logger = logging.getLogger(__name__)

    # Split on '/' to identify whether the image is bare, library-prefixed,
    # or registry-qualified.
    image_parts = image_name.split('/')

    if len(image_parts) == 1:
        # Bare form (e.g. "gcc:14") — official by convention.
        base_name = image_parts[0].split(':')[0]
        is_official = base_name in ALLOWED_OFFICIAL_IMAGES
    elif len(image_parts) == 2:
        # Either "library/gcc:14" (official) or "user/image:tag" (not).
        if image_parts[0] == 'library':
            base_name = image_parts[1].split(':')[0]
            is_official = base_name in ALLOWED_OFFICIAL_IMAGES
        else:
            is_official = False
    elif len(image_parts) == 3:
        # Fully-qualified form: only docker.io/library/<official> counts.
        if image_parts[0] == 'docker.io' and image_parts[1] == 'library':
            base_name = image_parts[2].split(':')[0]
            is_official = base_name in ALLOWED_OFFICIAL_IMAGES
        else:
            is_official = False
    else:
        is_official = False

    if not is_official:
        logger.warning(f"Image '{image_name}' is not in the allowed official images list")
        logger.warning(f"Allowed images: {', '.join(ALLOWED_OFFICIAL_IMAGES)}")

    return is_official


def pull_docker_image(image_name):
    """Run ``docker pull`` for ``image_name``; return True on success."""
    logger = logging.getLogger(__name__)

    logger.info(f"Pulling Docker image {image_name}...")
    try:
        # check=True raises on non-zero exit; pull progress streams to stdout.
        subprocess.run(["docker", "pull", image_name], check=True)
        logger.info(f"Successfully pulled {image_name}")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to pull image {image_name}")
        logger.error(f"Error: {e}")
        return False


def check_docker_image_exists(image_name):
    """Return True if ``image_name`` is available locally or in the registry."""
    # Local check first — avoids a network round-trip when the image is cached.
    cmd_local = ["docker", "image", "inspect", image_name]
    if subprocess.run(cmd_local, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode == 0:
        return True

    # Fall back to a remote manifest check.
    cmd_remote = ["docker", "manifest", "inspect", image_name]
    return subprocess.run(cmd_remote, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode == 0


def modify_dockerfile(dockerfile_path, config):
    """Substitute base image and ARGs in a Dockerfile template in-place.

    Constructs the ifsdata / climate / rtables URLs from the config, optionally
    HEAD-checks every data URL, then rewrites the Dockerfile placeholders.

    The replacements are tolerant of missing keys — if the CI Dockerfile drops
    e.g. ``OPENIFS_EXPT_URL`` and ``SCM_URL``, the corresponding ``.replace()``
    calls simply find no match and no-op. Callers can therefore share this
    function across the full and the CI Dockerfile templates.
    """
    logger = logging.getLogger(__name__)

    # Build-mode switch:
    # - True  (default): normal Docker image build includes OpenIFS runtime data
    # - False: CI/minimal image build skips ifsdata/climate/rtables URL handling
    include_openifs_data_downloads = config.get('include_openifs_data_downloads', True)

    openifs_version = config['openifs_version']
    climate_version = config.get('climate_version', '')

    if include_openifs_data_downloads:
        base_data_url = config['openifs_data_base_url']
        ifsdata_url = f"{base_data_url}/{openifs_version}/ifsdata/ifsdata.tar.gz"
        climate_url = f"{base_data_url}/{openifs_version}/{climate_version}/{openifs_version}_{climate_version}_159.tar.gz"
        rtables_url = f"{base_data_url}/{openifs_version}/rtables/rtables.tar.gz"

        logger.info("Constructed data URLs:")
        logger.info(f"  IFS data: {ifsdata_url}")
        logger.info(f"  Climate data: {climate_url}")
        logger.info(f"  RTables: {rtables_url}")
    else:
        ifsdata_url = ""
        climate_url = ""
        rtables_url = ""
        logger.info("Skipping OpenIFS data URL construction for CI/minimal Docker build")

    # Validate URLs before modifying the Dockerfile so a build doesn't run for
    # an hour and then fail on a wget. Only optional URLs (expt/SCM) may be
    # absent — those are skipped without failing.
    if not config.get('skip_url_validation', False):
        logger.info("Validating URLs to use in the Dockerfile, before image build...")

        urls_to_check = {
            'SCM package': config.get('scm_url', ''),
            'Experiment package': config.get('openifs_expt_url', ''),
        }

        if include_openifs_data_downloads:
            urls_to_check.update({
                'IFS data': ifsdata_url,
                'Climate data': climate_url,
                'RTables': rtables_url,
            })

        all_valid = True

        for name, url in urls_to_check.items():
            if not url:  # Optional URL not configured — skip.
                logger.warning(f"{name}: URL not configured, skipping check")
                continue

            accessible, error = check_url_accessible(url)
            if accessible:
                logger.info(f"{name}: {url} - Accessible")
            else:
                logger.error(f"{name}: {url} - Not accessible: {error}")
                all_valid = False

        if not all_valid:
            logger.error("Some URLs are not accessible - build will likely fail")
            logger.error("Set 'skip_url_validation: True' in config to bypass this check")
            sys.exit(1)

        logger.info("All data URLs are accessible")
    else:
        logger.warning("URL validation skipped (skip_url_validation=True)")
        logger.warning("Image build will fail if URL not available")

    # Rewrite the template. Each .replace() targets the placeholder exactly
    # as it appears in the Dockerfile (an empty ARG default or the FROM line).
    with open(dockerfile_path, "r") as file:
        content = file.read()

    content = content.replace('FROM docker.io/library/gcc:13.2.0-bookworm',
                              f'FROM docker.io/library/gcc:{config["base_docker_image"]}')
    content = content.replace('ARG OPENIFS_DIR=', f'ARG OPENIFS_DIR={config["openifs_version"]}')
    content = content.replace('ARG OPENIFS_EXPT_URL=', f'ARG OPENIFS_EXPT_URL={config.get("openifs_expt_url", "")}')
    content = content.replace('ARG SCM_URL=', f'ARG SCM_URL={config.get("scm_url", "")}')
    content = content.replace('ARG OPENIFS_REPO_URL=', f'ARG OPENIFS_REPO_URL={config["openifs_repo_url"]}')
    content = content.replace('ARG OPENIFS_BRANCH=', f'ARG OPENIFS_BRANCH={config["openifs_branch"]}')
    content = content.replace('ARG IFSDATA_URL=', f'ARG IFSDATA_URL={ifsdata_url}')
    content = content.replace('ARG CLIMATE_URL=', f'ARG CLIMATE_URL={climate_url}')
    content = content.replace('ARG CLIMATE_VERSION=', f'ARG CLIMATE_VERSION={climate_version}')
    content = content.replace('ARG RTABLES_URL=', f'ARG RTABLES_URL={rtables_url}')

    with open(dockerfile_path, 'w') as f:
        f.write(content)

    logger.info(f"Modified Dockerfile written to {dockerfile_path}")


def build_docker_image(dockerfile_path, image_name, build_dir, no_cache=False):
    """Run ``docker build`` for ``image_name`` against ``dockerfile_path``.

    The build context is ``build_dir`` (where the cloned OpenIFS source has
    been placed). ``no_cache=True`` is the safer default for CI — it forces a
    full rebuild, paying the cost in exchange for reproducibility.
    """
    logger = logging.getLogger(__name__)

    cmd = ["docker", "build"]
    if no_cache:
        cmd.append("--no-cache")
    cmd += ["-t", image_name, "-f", dockerfile_path, "."]

    logger.info(f"Executing image build using: {' '.join(cmd)}")
    subprocess.run(cmd, check=True, cwd=build_dir)
    logger.info("Docker image build completed")
