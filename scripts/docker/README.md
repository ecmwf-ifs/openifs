# OpenIFS Docker Builder

Automated Docker container creation for the stand-alone OpenIFS model.

## Quick Start

### Prerequisites

- Docker installed and running
- Python 3 with `gitpython` and `pyyaml` modules
- Git configured with SSH access to the OpenIFS repository

#### Python Environment Setup

It's recommended to use a virtual environment to install the required Python packages:

```bash
cd scripts/docker

# Create a virtual environment
python3 -m venv openifs-env

# Activate the virtual environment
# On macOS/Linux:
source openifs-env/bin/activate

# Install required packages
python3 -m pip install gitpython pyyaml

# Verify installation
python3 -c "import git, yaml; print('Packages installed successfully')"
```

**Note:** Keep the virtual environment activated when running the build script. To deactivate later, simply run `deactivate`.

### Basic Usage

1. Edit the configuration file:

   ```bash
   cp config/create_openifs_docker.yml config/my_config.yml
   ```

   Edit `my_config.yml` with your settings. In particular set `openifs_build_docker_dir` to an
   existing location.

1. Run the build script:

   ```bash
   python3 create-oifs-docker.py -c config/my_config.yml
   ```

The script will, depending on settings in the yaml config file:

- Clone the OpenIFS repository
- Copy SCM experiment data
- Build a Docker image with all dependencies
  - GCC compiler suite (version specified in config)
  - OpenMPI for parallel execution
  - NetCDF libraries for data I/O
  - LAPACK, Eigen3, and Boost libraries
  - Python 3 with required packages
  - OpenIFS source code and SCM test data
- Run tests to verify the installation

## Detailed Configuration

### Configuration File Options

The configuration file, `config/create_openifs_docker.yml`, controls the build of the image/container, clone of the OpenIFS repository and whether the tests are run or not. Below is more detail about the configuration options:

#### OpenIFS Settings

```yaml
# OpenIFS version (used for directory naming and image tagging)
openifs_version: "48r1"

# Git branch to extract from repository
openifs_branch: "main"

# Repository URL (requires SSH access)
openifs_repo_url: "git@github.com:ecmwf-ifs/openifs.git"

# SCM experiment data URL (tar.gz or tar file)
scm_url: https://openifs.ecmwf.int/data/scm/48r1/scm_openifs_48r1.tar.gz

# Clone repository (True) or use existing directory (False)
clone_openifs: True

# Force removal of existing clone without prompting
force_reclone: False

# Run openifs build command after building image
run_build: True

# Run openifs-test tests after building image
run_tests: True

# Run standard SCM tests - BOMEX, DYCOMS, TWPICE - after main tests
run_scm_test: False
```

#### Docker Settings

```yaml
# Base GCC Docker image version (e.g., "13", "13.2.0-bookworm")
base_docker_image: "13"

# Path to Dockerfile template
docker_template: "./Dockerfile"

# Force rebuild even if image exists
force_rebuild: True

# Force removal of test container after tests complete (True = remove with --rm, False = keep for inspection)
remove_test_container: False
```

#### Directory Settings

```yaml
# Directory for Docker build context and cloned repository
openifs_build_docker_dir: "~/oifs_docker_create_dir"
```

## Running the Container

The [Quick Start](#quick-start) section describes using `create-oifs-docker.py` to automate the installation, build and testing of OpenIFS in a container. If this is successful, the container is removed and success is reported, e.g.

```bash
[INFO] __main__.run_openifs_test : SUCCESS: All OpenIFS tests passed for openifs-48r1-gcc13:main
[INFO] __main__.main : All tests passed successfully
[INFO] __main__.main : ======================================================================
[INFO] __main__.main : Summary:
[INFO] __main__.main :   Image: openifs-48r1-gcc13:main
[INFO] __main__.main :   Built: Yes
[INFO] __main__.main :   Tests: Passed
[INFO] __main__.main : ======================================================================
```

The container is only removed if the scripts complete successfully. While the removal of the container saves space, all the build and test results are removed and irrecoverable.

Since the image has been built and is retained, the container can be re-created and started again using the following:

```bash
# Interactive shell
docker run -it <image-name>
```

where the image name can be taken from the report, e.g `openifs-48r1-gcc13:main`, or it can be found in the output of `docker ps -a` which lists stopped containers alongside their image names.

This is a clean container in which `source oifs-config.edit_me.sh` is run upon start up. This container can be used as an environment for building, testing and running OpenIFS. For example, once started `$OIFS_TEST/openifs-test.sh -cbt` can be executed (see [Main OpenIFS README](../../README.md)).

## `create-oifs-docker.py` workflow details

### Step 1: Validation

- Checks for required Python modules
- Validates base Docker image is from official sources (security)
- Checks if base image exists locally, pulls if needed

### Step 2: Repository Setup

- Shallow clones OpenIFS from specified branch (if `clone_openifs: True`)
- Copies SCM experiment data to build directory
- Updates configuration files with correct paths

### Step 3: Docker Build

- Creates Dockerfile from template with specified GCC version
- Installs all required dependencies
- Copies OpenIFS code and data into container
- Configures environment for non-root user

### Step 4: Testing (Optional)

- Runs OpenIFS test suite inside container
- Tests creation, compilation, and execution
- Reports success/failure with detailed logging

## Troubleshooting

### Image Already Exists

- Set `force_rebuild: True` to rebuild
- Or manually remove the image

### Clone Directory Exists

- Set `force_reclone: True` to remove and re-clone
- Or set `clone_openifs: False` to use existing directory

### Base Image Not Found

- Script will attempt to pull from Docker Hub
- Ensure Docker is running and you have internet access
- Verify the GCC version exists: https://hub.docker.com/_/gcc

### Test Failures

- Check log file in `docker_bld_logfiles/`
- Verify SCM data directory exists and contains required files
- Consider different GCC version if compatibility issues arise

## Security Notes

- Only official Docker images from approved sources are allowed
- Base images are validated before pulling
- Container runs as non-root user (uid 1000)

## Supported Configurations

Tested with:

- OpenIFS 48r1
- GCC 11.2, 12.2, 13.2
- Debian base images
