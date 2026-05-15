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


def _expand_strings(config):
    """Expand environment variables and ~ in all top-level string values."""
    for key, value in config.items():
        if isinstance(value, str):
            value = os.path.expandvars(value)   # expands $HOME, $VAR
            value = os.path.expanduser(value)   # expands ~
            config[key] = value
    return config


def main(config_path):
    """Read a YAML config file and return it as a dictionary.

    All top-level string values have ``$VAR`` and ``~`` expanded.

    Args
    ---
    config_path (str) : Path to the YAML configuration file.

    Return
    ------
    config (dict) : Configuration with all string values expanded.
    """
    import yaml

    config_path = os.path.abspath(config_path)

    with open(config_path, "r") as f:
        config = yaml.safe_load(f) or {}

    _expand_strings(config)

    return config


if __name__ == "__main__":
    main()
