#!/usr/bin/env python3

# (C) Copyright 2021- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

from datetime import datetime, timedelta
import argparse
import os

date = "20221106"
nsteps = 12
step = 1 # in hours
grid = "AV"
target = "ecwam_forcing.grib"

###############################################################################
# Overwrite defaults with command line arguments

parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""Generate MARS request to download ecwam forcings (10u, 10v, ci).

The forcing consists of
  - analysis at given <date>
  - 2 further analysis steps with 6 hours interval
  - forecast step 0 starting at last analysis step
  - <steps> further forecast steps with <step> hours interval

Note that the last analysis step and the first forecast step 0 are in fact identical.
The ecwam run is expected to ignore one of these

Example use:

      """+os.path.basename(__file__)+""" [OPTIONS] > request
      mars request""")

parser.add_argument("--date",   help="Date in format yyyymmdd [default="+date+"]", type=str)
parser.add_argument("--steps",  help="Number of steps [default="+str(nsteps)+"]", type=int)
parser.add_argument("--step",   help="Step in hours [default="+str(step)+"]", type=int)
parser.add_argument("--grid",   help="Output Grid [default=AV (archived value)]", type=str)
parser.add_argument("--target", help="Output grib file [default="+target+"]", type=str)
args = parser.parse_args()

if args.date:
  date = args.date
if args.steps:
  nsteps = args.steps
if args.step:
  step = args.step
if args.grid:
  grid = args.grid
if args.target:
  target = args.target

##############################################################################
# Datetime

def to_datetime(date):
  dtformat = "%Y%m%d%H%M%S"
  return datetime.strptime(date,dtformat[:len(date)-2])

def date_str(date):
  return date.strftime("%Y%m%d")

def time_str(dt):
  if isinstance(dt, datetime):
    return dt.strftime("%H")
  if isinstance(dt, timedelta):
    return str(int(dt.total_seconds())//3600)

def assemble_times(date, step, steps):
  times = {}
  for i in range(steps+1):
    time = ( date + i*step )
    if time.date() in times:
        times[time.date()].append(time)
    else:
        times[time.date()] = [time]
  return times

an_begin = to_datetime(date)
an_step  = timedelta(hours=6)
an_steps = 2
an_times  = assemble_times(date=an_begin, step=an_step, steps=an_steps)

fc_begin  = an_begin + an_steps * an_step
fc_step   = timedelta(hours=step)
fc_steps  = nsteps
fc_length = fc_steps * fc_step

##############################################################################
# Print mars request to stdout

def once():
  if not once.added:
    once.added = True
    return """
  PARAM=10U/10V/CI,
  STREAM=oper,
  CLASS=od,
  EXPVER=0001,
  REPRES=GG,
  TYPE=AN,
  LEVTYPE=SFC,
  LEVELIST=OFF,
  TARGET={target},
  GRID={grid},""".format(target=target,grid=grid)
  else:
    return ""
once.added = False

for date, datetimes in an_times.items():
  print("RETRIEVE,"+once()+"""
  DATE="""+date_str(date)+""",
  TIME="""+'/'.join([time_str(t) for t in datetimes]))

print("""RETRIEVE,
  TYPE=FC,
  DATE="""+date_str(fc_begin)+""",
  TIME="""+time_str(fc_begin)+""",
  STEP=0/to/"""+time_str(fc_length)+"/by/"+time_str(fc_step))
