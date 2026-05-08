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

import sys

def main(pymod_list) :

    """
    Check whether required Python modules are available for import.

    This function iterates through a list of expected Python modules and uses
    `importlib.util.find_spec()` to verify whether each can be imported without
    actually loading it. If any required module cannot be found, an error message
    is displayed and the program exits.

    Args:
        pymod_list (list[str]): A list of module names (as strings) to check.

    """

    import importlib.util

    module_avail = True

    for mod in pymod_list :
        if importlib.util.find_spec(mod) is not None:
            print(f"INFO : Module '{mod}' is available.")
        else:
            print(f"""
ERROR : {mod} is not available. Please check and re-run.
        Check README.md for info about how create and install
        in a virtual environment for specific modules. 
        """)
            module_avail = False

    if not module_avail :
        sys.exit()
    
if __name__ == "__main__":
    main()
