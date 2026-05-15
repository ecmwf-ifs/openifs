# scripts/shared — shared library

Generic utilities shared across the OpenIFS automation scripts:
`../bootstrap/docker/` and `../ci/docker_ci/`.

## Modules

| Module | Purpose |
| --- | --- |
| `shared_helpers.py` | Generic (non-Docker) helpers — `format_duration`, `timer`, `slug`, `shallow_clone`, `check_url_accessible`, `resolve_openifs_source`, `populate_from_local`, `patch_oifs_home`, `stage_branch_source` |
| `docker_lib.py` | Docker-specific helpers — base-image validate / pull / exists, Dockerfile ARG rewrite, image build |
| `ci_lib.py` | CI-specific helpers — report filename + control tarball name, synthetic-report writer, ctest-output appender, `build_test_commands`, `build_ci_summary` |
| `read_yml_config.py` | Load a YAML config file and expand `~` / `$VAR` in string values |
| `setup_logging.py` | Configure the standard `logging` module to write to both a log file and the console |
| `find_py_packages.py` | Check that a list of Python modules is importable; exits with a friendly error if any are missing |

## Importing from a driver script

Each driver inserts `scripts/shared/` onto `sys.path` at startup:

```python
_SHARED_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "shared")
if _SHARED_DIR not in sys.path:
    sys.path.insert(0, _SHARED_DIR)

import read_yml_config
import setup_logging
import find_py_packages
import shared_helpers
import docker_lib           # Docker drivers only
import ci_lib               # CI drivers only
```

`docker_lib` itself imports `check_url_accessible` from `shared_helpers`,
so the dependency direction is one-way: `docker_lib` -> `shared_helpers`,
never the reverse.  `ci_lib` imports from `shared_helpers` only.
