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
import os
import subprocess
import sys
import shutil
import git
import argparse
import logging
import time
from datetime import timedelta
from contextlib import contextmanager

import setup_logging
import read_yml_config
import find_py_packages

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

def format_duration(seconds):
    """Format duration in seconds to human-readable string."""
    return str(timedelta(seconds=int(seconds)))

@contextmanager
def timer(description, timings_dict, key):
    """Context manager for timing code blocks."""
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


def is_official_docker_image(image_name):
    """
    Check if an image is from an official/trusted source.
    
    Args:
        image_name: Image name (e.g., 'gcc:13.2.0-bookworm' or 'myuser/gcc:tag')
    
    Returns:
        bool: True if official, False otherwise
    """
    logger = logging.getLogger(__name__)
    
    # List of allowed official images (whitelist approach - most secure)
    ALLOWED_OFFICIAL_IMAGES = [
        'gcc',
        'ubuntu',
        'debian',
    ]
    
    # Extract image name without tag
    # Handle formats: gcc:tag, docker.io/library/gcc:tag, user/image:tag
    image_parts = image_name.split('/')
    
    if len(image_parts) == 1:
        # Format: gcc:tag (official image)
        base_name = image_parts[0].split(':')[0]
        is_official = base_name in ALLOWED_OFFICIAL_IMAGES
    elif len(image_parts) == 2:
        # Could be: library/gcc or user/image
        if image_parts[0] == 'library':
            base_name = image_parts[1].split(':')[0]
            is_official = base_name in ALLOWED_OFFICIAL_IMAGES
        else:
            # user/image format - not official
            is_official = False
    elif len(image_parts) == 3:
        # Format: docker.io/library/gcc
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
    """
    Pull a Docker image from registry.
    
    Args:
        image_name: Full image name with tag
    
    Returns:
        bool: True if successful, False otherwise
    """
    logger = logging.getLogger(__name__)
    
    logger.info(f"Pulling Docker image {image_name}...")
    pull_cmd = ["docker", "pull", image_name]
    
    try:
        # Show pull progress in real-time
        result = subprocess.run(pull_cmd, check=True)
        logger.info(f"Successfully pulled {image_name}")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to pull image {image_name}")
        logger.error(f"Error: {e}")
        return False

def check_docker_image_exists(image_name):
    """Check if Docker image exists locally or in registry."""
    logger = logging.getLogger(__name__)
    
    # Check locally first
    cmd_local = ["docker", "image", "inspect", image_name]
    result = subprocess.run(cmd_local, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if result.returncode == 0:
        return True
    
    # Check remote registry
    cmd_remote = ["docker", "manifest", "inspect", image_name]
    result = subprocess.run(cmd_remote, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return result.returncode == 0

def shallow_clone(repo_url, clone_dir, branch="main", force=False):
    """Shallow clones a repository at a specific branch, with default branch 'main'."""
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
    repo = git.Repo.clone_from(repo_url, clone_dir, depth=1, branch=branch)

    return

def check_url_accessible(url, timeout=10):
    """
    Check if a URL is accessible and returns a valid response.
    
    Args:
        url: URL to check
        timeout: Request timeout in seconds (default: 10)
    
    Returns:
        tuple: (accessible: bool, error_msg: str or None)
    """
    import urllib.request
    import urllib.error
    
    try:
        # Use HEAD request to check without downloading the file
        req = urllib.request.Request(url, method='HEAD')
        with urllib.request.urlopen(req, timeout=timeout) as response:
            if response.status == 200:
                return True, None
            else:
                return False, f"HTTP {response.status}"
    except urllib.error.HTTPError as e:
        return False, f"HTTP {e.code}: {e.reason}"
    except urllib.error.URLError as e:
        return False, f"URL Error: {e.reason}"
    except Exception as e:
        return False, f"Error: {str(e)}"

def modify_dockerfile(dockerfile_path, config):
    """
    Modify the Dockerfile to replace the base image version and set build arguments.
    Constructs and validates data file URLs based on config values.
    
    Args:
        dockerfile_path: Path to Dockerfile template
        config: Configuration dictionary containing all build parameters
    """
    logger = logging.getLogger(__name__)
    
    # Construct data file URLs
    openifs_version = config['openifs_version']
    climate_version = config['climate_version']
    base_data_url = config['openifs_data_base_url']
    
    ifsdata_url = f"{base_data_url}/{openifs_version}/ifsdata/ifsdata.tar.gz"
    climate_url = f"{base_data_url}/{openifs_version}/{climate_version}/{openifs_version}_{climate_version}_159.tar.gz"
    rtables_url = f"{base_data_url}/{openifs_version}/rtables/rtables.tar.gz"
    
    logger.info(f"Constructed data URLs:")
    logger.info(f"  IFS data: {ifsdata_url}")
    logger.info(f"  Climate data: {climate_url}")
    logger.info(f"  RTables: {rtables_url}")
    
    # Validate URLs before modifying Dockerfile
    if not config.get('skip_url_validation', False):
        logger.info("Validating URLs to use in the Dockerfile, before image build...")
        
        urls_to_check = {
            'IFS data': ifsdata_url,
            'Climate data': climate_url,
            'RTables': rtables_url,
            'SCM package': config.get('scm_url', ''),
            'Experiment package': config.get('openifs_expt_url', ''),
        }
        
        all_valid = True
        
        for name, url in urls_to_check.items():
            if not url:  # Skip empty URLs
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
    
    # Now modify the Dockerfile with validated URLs
    with open(dockerfile_path, "r") as file:
        content = file.read()

    # Replace placeholders
    content = content.replace('FROM docker.io/library/gcc:13.2.0-bookworm', 
                            f'FROM docker.io/library/gcc:{config["base_docker_image"]}')
    content = content.replace('ARG OPENIFS_DIR=', f'ARG OPENIFS_DIR={config["openifs_version"]}')
    content = content.replace('ARG OPENIFS_EXPT_URL=', f'ARG OPENIFS_EXPT_URL={config["openifs_expt_url"]}')
    content = content.replace('ARG SCM_URL=', f'ARG SCM_URL={config["scm_url"]}')
    content = content.replace('ARG OPENIFS_REPO_URL=', f'ARG OPENIFS_REPO_URL={config["openifs_repo_url"]}')
    content = content.replace('ARG OPENIFS_BRANCH=', f'ARG OPENIFS_BRANCH={config["openifs_branch"]}')
    content = content.replace('ARG IFSDATA_URL=', f'ARG IFSDATA_URL={ifsdata_url}')
    content = content.replace('ARG CLIMATE_URL=', f'ARG CLIMATE_URL={climate_url}')
    content = content.replace('ARG CLIMATE_VERSION=', f'ARG CLIMATE_VERSION={climate_version}')
    content = content.replace('ARG RTABLES_URL=', f'ARG RTABLES_URL={rtables_url}')
    
    with open(dockerfile_path, 'w') as f:
        f.write(content)

    logger.info(f"Modified Dockerfile written to {dockerfile_path}")
    
def update_oifs_home(oifs_config_path, openifs_version):
    """Update OIFS_HOME in oifs-config.edit_me.sh to use openifs_version.
    
    Returns:
        bool: True if file was modified or already correct, False if OIFS_HOME line not found.
    """
    from pathlib import Path
    logger = logging.getLogger(__name__)
    
    config_file = Path(oifs_config_path)
    if not config_file.exists():
        logger.error(f"Config file {oifs_config_path} not found")
        return False
    
    expected_line = f'export OIFS_HOME="${{HOME}}/{openifs_version}"'
    lines = config_file.read_text().splitlines()
    
    for i, line in enumerate(lines):
        if line.strip().startswith('export OIFS_HOME="${HOME}/'):
            if line.strip() == expected_line:
                logger.info(f"OIFS_HOME already correctly set in {oifs_config_path}")
                return True
            lines[i] = expected_line
            config_file.write_text('\n'.join(lines) + '\n')
            logger.info(f"Updated OIFS_HOME in {oifs_config_path}")
            return True
    
    logger.error(f"OIFS_HOME export line not found in {oifs_config_path}")
    return False

def build_docker_image(dockerfile_path, image_name, build_dir):
    """
    Builds a Docker image from the specified Dockerfile directory.
    By default, this is a clean build (includes no-cache), which is slower but safer
    """
    logger = logging.getLogger(__name__)
    
    cmd = ["docker", "build", "--no-cache", "-t", image_name, "-f", dockerfile_path, "."]

    logger.info(f"Executing image build using: {' '.join(cmd)}")

    subprocess.run(cmd, check=True, cwd=build_dir)

    logger.info("Docker image build completed")

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

    logger.info(f"Running openifs-test to create, build and test suite in container {image_name}...")
    logger.info("This may take 10-30 minutes depending on your system")

    # Build test command
    test_cmd = (
        f"source ~/{openifs_version}/oifs-config.edit_me.sh && "
        f"$OIFS_TEST/openifs-test.sh -cb -j 8''"
    )
    
    # Add OpenIFS tests to command if requested
    if run_tests :
        logger.info("OpenIFS test will be run after build success")
        test_cmd += " && $OIFS_TEST/openifs-test.sh -t"

    # Add SCM test if requested
    if run_scm_test:
        logger.info("SCM test will also be run after main tests")
        test_cmd += " && cd $OIFS_HOME && $SCM_TEST/callscm"
    
    # Build docker run command
    rm_flag = "--rm" if remove_container else ""    
    cmd = [
        "docker", "run", 
        *([rm_flag] if rm_flag else []),  # Add --rm only if specified, 
        image_name,
        "bash", "-lc",
        test_cmd
    ]

    logger.info(f"Running: {' '.join(cmd)}\n")
        
    try:
        result = subprocess.run(cmd)
        logger.info("OpenIFS built successfully")
        if run_tests:
            logger.info("OpenIFS tests passed successfully")
        if run_scm_test:
            logger.info("SCM test also passed successfully")
        if not remove_container:
            logger.info("Container was not removed. Use 'docker ps -a' to see it.")
        else : 
            logger.info("Container was removed")
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
        clone_repo = config.get('clone_openifs', True)

        if clone_repo:
            logger.info(f"Cloning OpenIFS repository to {openifs_dir}")
            shallow_clone(
                config['openifs_repo_url'], 
                openifs_dir, 
                branch=config['openifs_branch'],
                force=config.get('force_reclone', False)
            )
        else:
            if os.path.exists(openifs_dir):
                logger.info(f"Using existing OpenIFS repository at {openifs_dir}")
            else:
                logger.error(f"OpenIFS repository not found at {openifs_dir}")
                logger.error("Set 'clone_openifs: True' in config or provide existing directory")
                sys.exit(1)

    # Docker Image Build
    oifs_image_name = f"openifs-{config['openifs_version']}-gcc{config['base_docker_image']}:{config['openifs_branch']}"
    
    image_exists = check_docker_image_exists(oifs_image_name)
    force_rebuild = config.get('force_rebuild', False)
    
    should_build = False
    
    if not image_exists:
        logger.info(f"Docker image {oifs_image_name} does not exist - will build")
        should_build = True
    elif force_rebuild:
        logger.info(f"Docker image {oifs_image_name} exists but force_rebuild=True - will rebuild")
        should_build = True
    else:
        logger.info(f"Docker image {oifs_image_name} already exists - skipping build")
        logger.info("Set 'force_rebuild: True' in config to force rebuild")
    
    if should_build:
        with timer("Docker Image Build", timings, 'image_build'):
            logger.info(f"Building Docker image {oifs_image_name}...")
            build_docker_image(dockerfile_path, oifs_image_name, config['openifs_build_docker_dir'])
            logger.info(f"Docker image {oifs_image_name} built successfully!")
    else:
        timings['image_build'] = 0
    
    # OpenIFS Build and Test
    run_build = config.get('run_build', True)
    run_tests = config.get('run_tests', True)
    run_scm_test = config.get('run_scm_test', True)
    
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
                logger.error("Tests failed")
                if not should_build:
                    logger.error("Tests failed on existing image - consider setting 'force_rebuild: True'")
                else:
                    logger.error("Tests failed on newly built image - check build configuration")
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
    logger.info(f"  Built: {'Yes' if should_build else 'No (already exists)'}")
    logger.info(f"  OpenIFS Build: {'Passed' if run_build and test_success else 'Failed' if run_build else 'Skipped'}")
    logger.info(f"  OpenIFS Tests: {'Passed' if run_tests and test_success else 'Failed' if run_tests else 'Skipped'}")
    logger.info(f"  SCM Tests: {'Passed' if run_scm_test and test_success else 'Failed' if run_scm_test else 'Skipped'}")
    logger.info("=" * 70)
    logger.info("Timing Summary:")
    logger.info(f"  Image Validation:     {format_duration(timings['image_validation'])}")
    logger.info(f"  Dockerfile Prep:      {format_duration(timings['dockerfile_prep'])}")
    logger.info(f"  Repository Setup:     {format_duration(timings['repo_setup'])}")
    if should_build:
        logger.info(f"  Image Build:          {format_duration(timings['image_build'])}")
    else:
        logger.info(f"  Image Build:          Skipped")
    if run_build:
        logger.info(f"  Build & Test:         {format_duration(timings['build_and_test'])}")
    else:
        logger.info(f"  Build & Test:         Skipped")
    logger.info("  " + "-" * 66)
    logger.info(f"  Total:                {format_duration(total_time)}")
    logger.info("=" * 70)

if __name__ == "__main__":
    
    main()