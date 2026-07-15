#!/usr/bin/env python3

# (C) Copyright 2021- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

import argparse
import yaml
import os

parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""Read and extract information of parameters""")

parser.add_argument("parameters", help="List of parameter names", type=str, nargs='*' )
parser.add_argument("--request",  help="Return index for paramter", type=str, )
args = parser.parse_args()

parameters_file=os.path.dirname(os.path.realpath(__file__))+'/../parameters.yml'

f = open(parameters_file,'r')
yaml_document = f.read()
parameters = yaml.safe_load(yaml_document)['parameters']
f.close()

# print(parameters)


keys = args.parameters

if args.request=='index':
    index = [ parameters[key]['index'] if isinstance(key, str) else parameters[f'{key:03d}']['index'] for key in keys  ]
    print (' '.join([ str(idx) for idx in index]))
else:
    print('ERROR: unsupported value for --request')
    print('       Supported: --request=index')
