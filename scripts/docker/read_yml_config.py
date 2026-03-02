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

def main(config_path):

    """
    read_yaml_config uses pyyaml to read config_path and store the 
    contents in the dictionary config.

    Args 
    --- 
    config_path (str) : Path for the file yaml configuration file, e.g. $OIFS_DOCKER/config/create_openifs_config.yml.
  
    Return
    ------
    config (dictionary) : Contains the configuration that has been read from the the yaml file
    """

    # import yaml so config file can be read. 
    import yaml

    with open(config_path, "r") as f:
        config = yaml.safe_load(f)
    
    for key, value in config.items():
        if isinstance(value, str):
            value = os.path.expandvars(value)   # expands $HOME
            value = os.path.expanduser(value)   # expands ~ to home dir
            config[key] = value

    return config

if __name__ == "__main__":
    main()
