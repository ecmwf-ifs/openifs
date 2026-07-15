#!/usr/bin/env python3

# (C) Copyright 2021- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

import yaml
from datetime import datetime, timedelta
import sys
import re
import argparse

parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""Read and extract components from YAML configuration file""")

parser.add_argument("config",   help="YAML configuration file", type=str)
parser.add_argument("request",  help="request", type=str)

parser.add_argument("--format",  help="Format for key", type=str)
parser.add_argument("--default", help="Default value if key not found", type=str)

args = parser.parse_args()

f = open(args.config,'r')
yaml_document = f.read()
config = yaml.safe_load(yaml_document)
f.close()

request = args.request
default = args.default

def to_seconds(time):
    hours = 0
    minutes = 0
    seconds = 0
    try:
        if isinstance(time,str):
            hh_mm_ss = time.split(':')
            if len(hh_mm_ss) == 1:
                seconds = int(hh_mm_ss[0])
            if len(hh_mm_ss) == 2:
                hours   = int(hh_mm_ss[0])
                minutes = int(hh_mm_ss[1])
            if len(hh_mm_ss) == 3:
                hours   = int(hh_mm_ss[0])
                minutes = int(hh_mm_ss[1])
                seconds = int(hh_mm_ss[2])
        else:
                seconds = int(time)
    except ValueError as e:
        raise ValueError('Cannot convert "'+time+'" to seconds')
    return hours * 3600 + minutes * 60 + seconds

def to_datetime(date):
  dtformat = "%Y%m%d%H%M%S"
  return datetime.strptime(date,dtformat[:len(date)-2])

class Duration:
    def __init__(self,hours=0,minutes=0,seconds=0):
        self.total_seconds = hours * 3600 + minutes * 60 + seconds

    @classmethod
    def strptime(self,time,format):
        if format == "%H:%M:%S":
            match = re.match(r"([0-9]+):([0-9]+):([0-9]+)",time)
            if match:
                return Duration(hours=int(match.group(1)),minutes=int(match.group(2)),seconds=int(match.group(3)))
        elif format == "%H:%M":
            match = re.match(r"([0-9]+):([0-9]+)",time)
            if match:
                return Duration(hours=int(match.group(1)),minutes=int(match.group(2)))
        raise ValueError("Could not convert time "+str(time)+" to Duration with format "+format)
            

    def strftime(self,format):
        if "%H" in format:
            hours = self.total_seconds // 3600
            minutes = (self.total_seconds - hours * 3600) // 60
            seconds = self.total_seconds - hours * 3600 - minutes * 60
            return format.replace("%H",f'{hours:02d}').replace("%M",f'{minutes:02d}').replace("%S",f'{seconds:02d}')
        if "%M" in format:
            minutes = self.total_seconds // 60
            seconds = (self.total_seconds - minutes * 60) // 60
            return format.replace("%M",f'{minutes:02d}').replace("%S",f'{seconds:02d}')
        if "%S" in format:
            seconds = self.total_seconds
            return format.replace("%S",f'{seconds:02d}')

    def to_timedelta(self):
        return timedelta(seconds=self.total_seconds)

def flatten(l):
    flattened = []
    for item in l:
        if isinstance(item,list):
            for subitem in item:
                flattened.append(subitem)
        else:
            flattened.append(item)
    return flattened

def extract(config,keys,default=None):
    if not keys:
        return config

    key = keys[0]

    if isinstance(config,list) and key == 'size':
        return len(config)

    match_array = re.match(r"(.+)\[(.*)\]$",key)

    if match_array:
        key = match_array.group(1)
        access = match_array.group(2)
        if not key in config:
            if default is None:
                raise KeyError('Key '+key+' not found in '+str(config))
            array = default.split()
        else:
            array = config[key]
        if access == ":":
            return flatten([ extract(nested,keys[1:],default) for nested in array ])
        if re.match("([0-9]+)?-?[0-9]+",access):
            idx = int(eval(access))
            if default is not None and idx >= len(array):
                return default
            return extract(array[idx], keys[1:], default)
        raise ValueError("Could not match "+match_array.group(0)+ " as supported array syntax")

    if key in config:
        if len(keys) == 1:
            value = config[key]
            match_indirection = re.match(r"^\$\{(.*)\}$",str(value))
            if match_indirection:
                return extract(root,[match_indirection.group(1)],default)
            return value
        return extract(config[key], keys[1:], default)
    else:
        if key == 'time' and 'begin' in config and 'end' in config and 'timestep' in config:
            begin = config['begin']
            end = config['end']
            timestep = None
            for format in ("%H:%M","%H:%M:%S"):
                try:
                    timestep = Duration.strptime(config['timestep'],format).to_timedelta()
                except ValueError: pass
            if not timestep:
                raise ValueError('Could not parse timestep '+config['timestep'])
            times = []
            time = begin
            while time < end:
                times.append(time)
                time += timestep
            return times

        if default is None:
            raise KeyError('Key '+key+' not found in '+str(config))
        return default

root = config
extracted = extract(config, request.split('.'), default)

def format(value, format):
        if not format:
            return str(value)
        if "%Y" in format:
            if isinstance(value,str):
                try:
                    value = to_datetime(value)
                except ValueError: pass
            if not isinstance(value, datetime):
                raise ValueError('Value is not a valid datetime: '+value)
            return value.strftime(format)
        if "%H" in format or "%M" in format or "%S" in format:
            return Duration(seconds=to_seconds(value)).strftime(format)
        if format == "seconds":
            return str(to_seconds(value))
        if format == "hours":
            return str(to_seconds(value)//3600)

        raise ValueError('Unrecognized format "'+format+'"')


### Format output
try:
    if isinstance(extracted, list):
        print( ' '.join([ format(value,args.format) for value in extracted ]))
    else:
        print(format(extracted,args.format))
except ValueError as e:
    if default is None:
        print("ERROR:",e)
    else:
        if isinstance(extracted, list):
            array = default.split()
            print( ' '.join([ value for value in array ]))
        else:
            print(default)
