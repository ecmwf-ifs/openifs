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
import re

parser = argparse.ArgumentParser()
parser.add_argument("grid_name", help="Named grid, e.g. O1280")
parser.add_argument("--grid_description", help="Output grid_description as input to create_wam_bathymetry and preproc", action='store_true')
parser.add_argument("--variable", help="Output variable")
args = parser.parse_args()

if not args.grid_name:
    print("ERROR: No grid_name argument provided")
    exit(1)

grid_name = args.grid_name

result = re.match("^O([0-9]+)$",grid_name)
if not result:
    print("ERROR: Could not match grid to octahedral grid name")
    exit(1)

N = int(result.group(1))

# Dictionary or first latitude of Gaussian grid,
# computed via:
#
#    for N in 16 24 32 48 64 80 96 128 160 200 256 320 400 512 576 640 800 1024 1280 1600 2000 4000 8000
#    do 
#        echo "lat0[$N]=$(atlas-gaussian-latitudes -N ${N} | head -1 | awk '{ print $2 }')"
#    done

lat0 = {}
lat0[16]=85.760587120444
lat0[24]=87.159094555863
lat0[32]=87.863798839233
lat0[48]=88.572168514007
lat0[64]=88.927735352296
lat0[80]=89.141519426461
lat0[96]=89.284227532514
lat0[128]=89.462821568577
lat0[160]=89.570089550607
lat0[200]=89.655964246870
lat0[256]=89.731148618413
lat0[320]=89.784876907219
lat0[400]=89.827874645894
lat0[512]=89.865508687700
lat0[576]=89.880445682778
lat0[640]=89.892396445590
lat0[800]=89.913910432567
lat0[1024]=89.932737928460
lat0[1280]=89.946187715666
lat0[1600]=89.956948491058
lat0[2000]=89.965557716640
lat0[4000]=89.982777782041
lat0[8000]=89.991388621915

if not N in lat0.keys():
    print("ERROR: Could not find first latitude in dictionary of Gaussian latitudes")
    exit(1)

ny = 2*N
dlon=360./(20+(N-1)*4)
north = lat0[N]
south = -lat0[N]
west = 0
east = 360.-dlon
iper = 1
irgg = 1

if args.variable:
    var=args.variable
    if var == "resol":
        print(N)
    if var == "cwamresol":
        print("{:05d}".format(N))
    if var == "north":
        print(north)
    if var == "south":
        print(south)
    if var == "west":
        print(west)
    if var == "east":
        print(east)
    if var == "dlat":
        print((north-south)/(ny-1))
    if var == "dlon":
        print(360./(20+(N-1)*4))
    if var == "irgg":
        print(1)
    exit(0)

if args.grid_description:
    print(N)                                  # resol
    print(north)                              # amonop
    print(south)                              # amosop
    print(west)                               # amowep
    print(east)                               # amoeap
    print(iper)                               # iper
    print(irgg)                               # irgg
    print(ny)                                 # ny
    for j in range(N):                        # pl
        print(20+j*4)
    for j in range(N):
        print(20+(N-1)*4 - j*4)

# Format of grid_description, copied from ifs-scripts/wave/wave_create_bathymetry:
#
#   Nj=$(grib_get -p Nj $file)
#   cwamresol=$(((Nj-1)))
#   echo $(((Nj-1))) > grid_description
#   echo $(grib_get -p latitudeOfFirstGridPointInDegrees $file) >> grid_description
#   echo $(grib_get -p latitudeOfLastGridPointInDegrees  $file) >> grid_description
#   echo $(grib_get -p longitudeOfFirstGridPointInDegrees $file) >> grid_description
#   echo $(grib_get -p longitudeOfLastGridPointInDegrees $file) >> grid_description
#   echo $(grib_get -p global $file) >> grid_description
#   echo $(grib_get -p global $file) >> grid_description
#   echo $Nj >> grid_description
#   # get the number of points per latitude
#   cat > filter <<EOF
#   print "[pl' ']" ;
#   EOF
#   pl=$(grib_filter filter $file)
#   for number in $pl
#   do
#     echo $number >> grid_description
#   done

