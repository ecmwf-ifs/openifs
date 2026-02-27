import os
import subprocess
import sys
import shutil
import git
import argparse
import logging

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

def modify_dockerfile(dockerfile_path, base_docker_version, openifs_dir, repo_url, branch, scm_url):
    """Modify the Dockerfile to replace the base image version and set OPENIFS_DIR."""
    logger = logging.getLogger(__name__)
    
    with open(dockerfile_path, "r") as file:
        content = file.read()

    # Replace placeholders
    content = content.replace('FROM docker.io/library/gcc:13.2.0-bookworm', 
                            f'FROM docker.io/library/gcc:{base_docker_version}')
    content = content.replace('ARG OPENIFS_DIR=', f'ARG OPENIFS_DIR={openifs_dir}')
    content = content.replace('ARG SCM_URL=', f'ARG SCM_URL={scm_url}')
    content = content.replace('ARG OPENIFS_REPO_URL=', f'ARG OPENIFS_REPO_URL={repo_url}')
    content = content.replace('ARG OPENIFS_BRANCH=', f'ARG OPENIFS_BRANCH={branch}')
    
    with open(dockerfile_path, 'w') as f:
        f.write(content)
    
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

def run_openifs_test(openifs_version, image_name):
    """Run openifs-test inside the Docker container and report results."""
    logger = logging.getLogger(__name__)

    logger.info(f"Running openifs-test to create, build and test suite in container {image_name}...")
    logger.info("This may take 10-30 minutes depending on your system")
        
    cmd = [
        "docker", "run", "--rm", image_name,
        "bash", "-lc",
        f"source ~/{openifs_version}/oifs-config.edit_me.sh && "
        f"$OIFS_TEST/openifs-test.sh -cbt -j 8 --arch=''"
    ]

    logger.info(f"Running: {' '.join(cmd)}\n")
        
    try:
        result = subprocess.run(cmd)
    except KeyboardInterrupt:
        logger.warning("Interrupted by user with ctrl-c")
        return False
        
    logger.info("")  # Add blank line after output
    
    if result.returncode != 0:
        logger.error(f"Build and/or test FAILED (exit code: {result.returncode})")
        return False
    else :
        # All steps passed
        logger.info(f"SUCCESS: All OpenIFS tests passed for {image_name}")
        return True
    
def main():
    
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

    # Check if base Docker image exists before proceeding
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

    docker_file_name = f"Dockerfile_{config['openifs_version']}_{config['base_docker_image']}"

    # Copy and Modify the Dockerfile
    dockerfile_path = os.path.join(config['openifs_build_docker_dir'], docker_file_name)

    # Check if Dockerfile exists and create backup
    if os.path.exists(dockerfile_path):
        logger.warning(f"Dockerfile {dockerfile_path} already exists, creating backup")
        shutil.copyfile(dockerfile_path, f"{dockerfile_path}.bak")
    else:
        logger.info(f"Creating Dockerfile {dockerfile_path}")

    # copy Dockefile template to dockerfile_path
    shutil.copyfile(config['docker_template'], dockerfile_path)
    modify_dockerfile(dockerfile_path, 
                      config['base_docker_image'], 
                      config['openifs_version'], 
                      config['openifs_repo_url'],
                      config['openifs_branch'],
                      config['scm_url']
                      )

    
    # Check for an existing OpenIFS directory in the docker build directory
    openifs_dir = os.path.join(config['openifs_build_docker_dir'], config['openifs_version'])

    if config['clone_openifs']:
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

    # Set the docker image name
    oifs_image_name = f"openifs-{config['openifs_version']}-gcc{config['base_docker_image']}:{config['openifs_branch']}"

    # Determine if Docker image needs to be built
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
    
    #Build Docker image if needed
    if should_build:
        logger.info(f"Building Docker image {oifs_image_name}...")
        build_docker_image(dockerfile_path, oifs_image_name, config['openifs_build_docker_dir'])
        logger.info(f"Docker image {oifs_image_name} built successfully!")
    
    # Run OpenIFS tests in docker container if configured
    run_tests = config.get('run_tests', True)
    
    if run_tests:
        logger.info("=" * 70)
        logger.info("Running openifs-test -cbt in container...")
        logger.info("=" * 70)
        
        test_success = run_openifs_test(config['openifs_version'], oifs_image_name)
        
        if test_success:
            logger.info("All tests passed successfully")
        else:
            logger.error("Tests failed")
            if not should_build:
                logger.error("Tests failed on existing image - consider setting 'force_rebuild: True'")
            else:
                logger.error("Tests failed on newly built image - check build configuration")
    else:
        logger.info("Skipping tests (run_tests: False in config)")
    
    # ===================================================================
    # SECTION 4: Summary
    # ===================================================================
    
    logger.info("=" * 70)
    logger.info("Summary:")
    logger.info(f"  Image: {oifs_image_name}")
    logger.info(f"  Built: {'Yes' if should_build else 'No (already exists)'}")
    logger.info(f"  Tests: {'Passed' if run_tests and test_success else 'Failed' if run_tests else 'Skipped'}")
    logger.info("=" * 70)

if __name__ == "__main__":
    
    main()