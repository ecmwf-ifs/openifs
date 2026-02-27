# OpenIFS Docker Builder

Automated Docker container creation for the stand-alone OpenIFS model.

## Quick Start

### Prerequisites

- Docker installed and running
- Python 3 with `git` and `yaml` modules
- Git configured with SSH access to the OpenIFS repository

#### Python Environment Setup

It's recommended to use a virtual environment to install the required Python packages:

```bash
# Create a virtual environment
python3 -m venv openifs-env

# Activate the virtual environment
# On macOS/Linux:
source openifs-env/bin/activate
# On Windows:
# openifs-env\Scripts\activate

# Install required packages
pip install gitpython pyyaml

# Verify installation
python3 -c "import git, yaml; print('Packages installed successfully')"
```

**Note:** Keep the virtual environment activated when running the build script. To deactivate later, simply run `deactivate`.

### Basic Usage

1. Clone this repository:

```bash
git clone <repository-url>
cd openifs-docker
```

1. Edit the configuration file:

```bash
cp config/create_openifs_docker.yml config/my_config.yml
# Edit my_config.yml with your settings
```

1. Run the build script:

```bash
python3 create-oifs-docker.py -c config/my_config.yml
```

The script will:

- Clone the OpenIFS repository
- Copy SCM experiment data
- Build a Docker image with all dependencies
- Run tests to verify the installation

## Detailed Configuration

### Configuration File Options

Edit `config/create_openifs_docker.yml` to customize your build:

#### OpenIFS Settings

```yaml
# OpenIFS version (used for directory naming and image tagging)
openifs_version: "48r1"

# Git branch to extract from repository
openifs_branch: "main"

# Repository URL (requires SSH access)
openifs_repo_url: "git@github.com:ecmwf-ifs/openifs.git"

# Clone repository (True) or use existing directory (False)
clone_openifs: True

# Force removal of existing clone without prompting
force_reclone: False
```

#### Docker Settings

```yaml
# Base GCC Docker image version (e.g., "13", "13.2.0-bookworm")
base_docker_image: "13"

# Path to Dockerfile template
docker_template: "./Dockerfile"

# Force rebuild even if image exists
force_rebuild: True
```

#### Directory Settings

```yaml
# Directory for Docker build context and cloned repository
openifs_build_docker_dir: "~/oifs_docker_create_dir"

# SCM experiment data directory (must exist on your system)
scm_exp_datadir: "~/openifs-expt/scm_openifs"
```

#### Testing

```yaml
# Run tests after building image (creates, builds, and tests OpenIFS)
run_tests: True
```

### Command Line Options

```bash
python3 create-oifs-docker.py --config <yaml_file>
# or
python3 create-oifs-docker.py -c <yaml_file>
```

## What Gets Built

The Docker image includes:

- GCC compiler suite (version specified in config)
- OpenMPI for parallel execution
- NetCDF libraries for data I/O
- LAPACK, Eigen3, and Boost libraries
- Python 3 with required packages
- OpenIFS source code and SCM test data

## Output and Logs

Build logs are saved to:

```bash
<openifs_build_docker_dir>/docker_bld_logfiles/log_<version>_<gcc_version>.log
```

## Running the Container

After successful build, run the container:

```bash
# Interactive shell
docker run -it 48r1-gcc13:main bash

# Inside container, OpenIFS environment is pre-configured
source ~/48r1/oifs-config.edit_me.sh
cd $OIFS_EXPT
```

## Workflow Details

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

## License

[Add appropriate license information]

## Support

For issues and questions:

- Check log files for detailed error messages
- Verify all prerequisites are installed
- Ensure repository access is configured correctly