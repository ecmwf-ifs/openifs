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
import argparse
import logging
import os
import shutil
import subprocess
import sys
import time

# Generic helpers and Docker helpers all live in scripts/shared/.
_SHARED_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "shared")
if _SHARED_DIR not in sys.path:
    sys.path.insert(0, _SHARED_DIR)

import setup_logging
import read_yml_config
import find_py_packages

# Docker-specific helpers from scripts/shared/docker_lib.py.
from docker_lib import (  # type: ignore[import-not-found]
    build_docker_image,
    check_docker_image_exists,
    is_official_docker_image,
    modify_dockerfile,
    pull_docker_image,
)
# Generic (non-Docker) helpers from scripts/shared/shared_helpers.py.
from shared_helpers import (  # type: ignore[import-not-found]
    format_duration,
    move_to_backup,
    resolve_openifs_source,
    shallow_clone,
    slug,
    timer,
)

def parse_arguments() :
    parser = argparse.ArgumentParser(
        description=f"""
create_openifs_docker and the associated modules creates a 
container for the stand-alone package for OpenIFS. 

This script automates:
  1. Cloning OpenIFS from the specified branch
  2. Copying SCM experiment data
  3. Building a Docker image with GCC and required libraries
  4. Running OpenIFS tests to verify the installation

For detailed documentation, see README.md

Prerequisites:
  - Docker installed and running
  - Python 3 with git, yaml modules (see README.md for setup)
  - SSH access to OpenIFS repository

Usage:
    python3 create-oifs-docker.py -c config/create_openifs_docker.yml

For more information: README.md#detailed-configuration

""", 
       formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument("--config", "-c", type=str, 
                        help="YAML configuration file (see config/create_openifs_docker.yml)")
 
    args = parser.parse_args()  

    ######### Check for command line arguments ###########################################
    #
    # Check that user has provided a branch name, if not exit
    #  
    if args.config is None :
        parser.print_help()
        print(f"""
[ERROR]: User must provide an a yml config file using --config, e.g.
        <path_to_script>/create_openifs_driver.py -c config/create_openifs_config.yml
        """)
        sys.exit()
    
    ########################################################################################

    return args


def run_openifs_test(openifs_version, image_name,
                     run_tests=True, 
                     run_scm_test=True, 
                     remove_container=True):
    """
    Run openifs-test build inside the Docker container and report results.
    Tests are also run, depending on the arguments and the yml config
    
    Args:
        openifs_version: OpenIFS version string
        image_name: Docker image name to test
        run_tests : Run the OpenIFS tests
        run_scm_test : Run the standard SCM cases
        remove_container: If True, remove container after test completes (default: True)
    """
    logger = logging.getLogger(__name__)

    container_name = f"oifs-{openifs_version}"

    # Remove any existing container with the same name
    check_result = subprocess.run(
        ["docker", "inspect", "--format", "{{.Name}}", container_name],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
    )
    if check_result.returncode == 0:
        logger.warning(f"Container '{container_name}' already exists and will be removed")
        subprocess.run(["docker", "rm", "-f", container_name], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        logger.info(f"Existing container '{container_name}' removed successfully")
    
    # Start container with /bin/bash as main process
    logger.info(f"Starting container '{container_name}' from image {image_name}...")
    run_cmd = [
        "docker", "run", "-dit",
        "--name", container_name,
        image_name,
        "/bin/bash"
    ]
    subprocess.run(run_cmd, check=True)
    logger.info(f"Container '{container_name}' started. Re-enter later with:")
    logger.info(f"  docker start {container_name} && docker exec -it {container_name} /bin/bash")

    # Build test command (unchanged)
    test_cmd = (
        f"source ~/{openifs_version}/oifs-config.edit_me.sh && "
        f"$OIFS_TEST/openifs-test.sh -cb -j 8"
    )
    if run_tests:
        test_cmd += " && $OIFS_TEST/openifs-test.sh -t"
    if run_scm_test:
        test_cmd += " && cd $OIFS_HOME && $SCM_TEST/callscm"

    # Execute test command inside the running container via exec
    exec_cmd = [
        "docker", "exec", "-it",
        container_name,
        "bash", "-lc",
        test_cmd
    ]

    logger.info(f"Running tests via exec: {' '.join(exec_cmd)}\n")

    try:
        subprocess.run(exec_cmd, check=True)
        logger.info("OpenIFS built successfully")
        if run_tests:
            logger.info("OpenIFS tests passed successfully")
        if run_scm_test:
            logger.info("SCM test also passed successfully")
        if remove_container:
            subprocess.run(["docker", "rm", "-f", container_name], check=True)
            logger.info(f"Container '{container_name}' removed")
        else:
            logger.info(f"Container '{container_name}' left running. Use 'docker ps' to see it.")
            logger.info(f"Container can be restarted using 'docker exec -it {container_name} /bin/bash'")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"OpenIFS tests failed: {e}")
        logger.error(f"stdout: {e.stdout}")
        logger.error(f"stderr: {e.stderr}")
        if not remove_container:
            logger.info("Container was not removed. Use 'docker ps -a' to inspect it.")
        return False

def main():
    
    script_start_time = time.time()
    timings = {}
    
    # Read yaml config path from the command line
    cli_args = parse_arguments()

    # As the command line arguments have been accepted, now 
    # check that the "non-standard" python modules are available
    pymod_list=["git","yaml"]
    #
    find_py_packages.main(pymod_list)

    config = read_yml_config.main(cli_args.config)

    # Resolve openifs_source once and reuse the result everywhere.
    # The raw value comes from openifs_source (preferred) or the legacy
    # openifs_branch key.  resolve_openifs_source maps "" -> auto-detect.
    _raw = config.get('openifs_source') or config.get('openifs_branch', '')
    _source_kind, _source_value = resolve_openifs_source(_raw, __file__)
    if _source_kind == 'remote':
        _tag = slug(_source_value)
    else:
        _tag = f"local-{slug(os.path.basename(os.path.realpath(_source_value)))}"
    config.setdefault('openifs_branch', _tag)

    log_dir = os.path.join(config['openifs_build_docker_dir'], "docker_bld_logfiles")
    
    # Create directory if it doesn't exist
    os.makedirs(log_dir, exist_ok=True)

    log_file_path = os.path.join(log_dir, f"log_{config['openifs_version']}_{config['base_docker_image']}.log")

    # Setup to write logfile in the current working directory. Using default log info
    setup_logging.main(log_file_path)
    logger = logging.getLogger(__name__)

    # Docker Base Image Validation
    with timer("Docker Base Image Validation", timings, 'image_validation'):
        base_image = f"gcc:{config['base_docker_image']}"
        
        # Security check: only allow official/vetted images
        logger.info(f"Validating base Docker image {base_image}...")
        if not is_official_docker_image(base_image):
            logger.error(f"Security check failed: '{base_image}' is not an approved official image")
            logger.error("Only official Docker images are allowed for security reasons")
            logger.error("If you need to use a different image, add it to ALLOWED_OFFICIAL_IMAGES in the code")
            sys.exit(1)
        
        logger.info(f"Security check passed: {base_image} is an official image")
        
        # Check if image exists locally
        logger.info(f"Checking if base Docker image {base_image} exists locally...")
        if not check_docker_image_exists(base_image):
            logger.warning(f"Base Docker image {base_image} not found locally")
            logger.info("Attempting to pull from Docker Hub...")
            
            if not pull_docker_image(base_image):
                logger.error(f"Failed to pull base Docker image {base_image}")
                logger.error("Please check your internet connection and Docker Hub status")
                logger.error(f"You can try manually: docker pull {base_image}")
                sys.exit(1)
        else:
            logger.info(f"Base Docker image {base_image} is available locally")

    # Dockerfile Preparation
    with timer("Dockerfile Preparation", timings, 'dockerfile_prep'):
        docker_file_name = f"Dockerfile_{config['openifs_version']}_{config['base_docker_image']}"
        dockerfile_path = os.path.join(config['openifs_build_docker_dir'], docker_file_name)

        # Check if Dockerfile exists and create backup
        if os.path.exists(dockerfile_path):
            logger.warning(f"Dockerfile {dockerfile_path} already exists, creating backup")
            shutil.copyfile(dockerfile_path, f"{dockerfile_path}.bak")
        else:
            logger.info(f"Creating Dockerfile {dockerfile_path}")

        # Check if template exists
        docker_template = config['docker_template']
        if not os.path.exists(docker_template):
            logger.error(f"Docker template file not found: {docker_template}")
            logger.error("Please check 'docker_template' path in your config file")
            sys.exit(1)

        shutil.copyfile(docker_template, dockerfile_path)
        modify_dockerfile(dockerfile_path, config)

    # OpenIFS Repository Setup
    with timer("OpenIFS Repository Setup", timings, 'repo_setup'):
        openifs_dir = os.path.join(config['openifs_build_docker_dir'], config['openifs_version'])

        # Reuse the resolution computed at the top of main() so we don't
        # accidentally re-interpret the derived openifs_branch tag as a
        # remote branch name (the bug that caused clone of
        # "local-openifs-casim" when openifs_source was empty).
        source_kind, source_value = _source_kind, _source_value

        if source_kind == 'remote':
            logger.info(f"Cloning branch '{source_value}' to {openifs_dir}")
            shallow_clone(
                config['openifs_repo_url'],
                openifs_dir,
                branch=source_value,
                force=config.get('force_reclone', False),
            )
            source_tag = slug(source_value)
        else:
            # 'auto' or 'local' — stage the resolved local path into the build dir.
            local_src = source_value
            if not os.path.isdir(local_src):
                logger.error(f"OpenIFS source not found at {local_src}")
                logger.error("Check 'openifs_source' in your config or re-run from inside the checkout")
                sys.exit(1)
            if os.path.abspath(local_src) == os.path.abspath(openifs_dir):
                logger.info(f"Source and build dir are the same ({openifs_dir}); skipping copy")
            elif os.path.exists(openifs_dir) and not config.get('force_reclone', False):
                logger.info(f"Using existing staged source at {openifs_dir} (force_reclone=False)")
            else:
                if os.path.exists(openifs_dir):
                    # force_reclone=True path: rename the existing staged
                    # tree to a timestamped backup so any uncommitted or
                    # unpushed work survives instead of being deleted.
                    move_to_backup(openifs_dir)
                logger.info(f"Copying local source {local_src} -> {openifs_dir}")
                shutil.copytree(
                    local_src, openifs_dir,
                    symlinks=True,
                    ignore=shutil.ignore_patterns(
                        '.git', 'build', '.cache', '.bootstrap',
                        '__pycache__', '*.pyc', '*.pyo', 'openifs-env',
                    ),
                )
            source_tag = f"local-{slug(os.path.basename(os.path.realpath(local_src)))}"

        # Propagate source_tag so modify_dockerfile and the image name both
        # reflect where this build came from (mirrors ci-oifs-docker.py).
        config['openifs_branch'] = source_tag

    # Docker Image Build
    oifs_image_name = f"openifs-{config['openifs_version']}-gcc{config['base_docker_image']}:{config['openifs_branch']}"
    
    force_rebuild = config.get('force_rebuild', False)

    with timer("Docker Image Build", timings, 'image_build'):
        logger.info(f"Building Docker image {oifs_image_name}...")
        if force_rebuild:
            logger.info("force_rebuild=True: building without cache")
        else:
            logger.info("force_rebuild=False: building with cache")
        logger.info(f"Building Docker image {oifs_image_name}...")
        build_docker_image(dockerfile_path, oifs_image_name, config['openifs_build_docker_dir'], no_cache=force_rebuild)
        logger.info(f"Docker image {oifs_image_name} built successfully!")
    
    # OpenIFS Build and Test
    run_build = config.get('run_build', True)
    run_tests = config.get('run_tests', True)
    run_scm_test = config.get('run_scm_test', True)
    test_success = False

    if run_build:
        with timer("OpenIFS Build and Test", timings, 'build_and_test'):
            test_success = run_openifs_test(
                config['openifs_version'],
                oifs_image_name,
                run_tests,
                run_scm_test,
                config.get('remove_test_container', True),
            )

            if test_success:
                logger.info("All tests passed successfully")
            else:
                logger.error("Tests failed - check build configuration")
    else:
        logger.info("Skipping build and tests (run_build: False in config)")
        timings['build_and_test'] = 0
    
    # Final Summary
    total_time = time.time() - script_start_time
    
    logger.info("=" * 70)
    logger.info("FINAL SUMMARY")
    logger.info("=" * 70)
    logger.info("Configuration:")
    logger.info(f"  Image: {oifs_image_name}")
    logger.info(f"  Cache: {'disabled (--no-cache)' if force_rebuild else 'enabled'}")
    logger.info(f"  OpenIFS Build: {'Passed' if run_build and test_success else 'Failed' if run_build else 'Skipped'}")
    logger.info(f"  OpenIFS Tests: {'Passed' if run_tests and test_success else 'Failed' if run_tests else 'Skipped'}")
    logger.info(f"  SCM Tests: {'Passed' if run_scm_test and test_success else 'Failed' if run_scm_test else 'Skipped'}")
    logger.info("=" * 70)
    logger.info("Timing Summary:")
    logger.info(f"  Image Validation:     {format_duration(timings['image_validation'])}")
    logger.info(f"  Dockerfile Prep:      {format_duration(timings['dockerfile_prep'])}")
    logger.info(f"  Repository Setup:     {format_duration(timings['repo_setup'])}")
    logger.info(f"  Image Build:          {format_duration(timings['image_build'])}")

    if run_build:
        logger.info(f"  Build & Test:         {format_duration(timings['build_and_test'])}")
    else:
        logger.info(f"  Build & Test:         Skipped")
    logger.info("  " + "-" * 66)
    logger.info(f"  Total:                {format_duration(total_time)}")
    logger.info("=" * 70)

if __name__ == "__main__":
    
    main()