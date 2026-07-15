# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

###############################################################################

# find FORD Fortran documentation generator
find_program( FORD_EXECUTABLE ford QUIET DOC "Fortran documentation generator" )

if( FORD_EXECUTABLE )
  set( FORD_FOUND TRUE )
endif()

mark_as_advanced(FORD_EXECUTABLE)
