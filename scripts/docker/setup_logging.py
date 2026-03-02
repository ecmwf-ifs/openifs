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

import logging

def main( logfile, level=logging.INFO ):

    """
    This function sets up a logging system that writes messages to both 
    a log file and the console. It defines a uniform message format that 
    includes the log level, logger name, function name, and message text.

    Args 
    --- 
    logfile (str) : Path to the file where log messages will be written. The file will be overwritten each time the function is called.
    level (string or integer), optional : default is INFO, which is quite verbose level options are INFO, WARNING and ERROR (in this implementation)
    """

    logging.basicConfig(
        level=level,  # Set the default log level to INFO
        format='[%(levelname)s] %(name)s.%(funcName)s : %(message)s',
        # set up 2 handlers so that output written to file and screen
        handlers=[
            logging.FileHandler(logfile, mode='w', encoding='utf-8'),  # Log to a file
            logging.StreamHandler()        # Log to screen
        ]
    )

if __name__ == "__main__":
    main()

