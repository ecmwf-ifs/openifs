# (C) Copyright 2022- ECMWF.
# (C) Copyright 2022- Meteo-France.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

include(FetchContent)

FetchContent_Declare(
   fypp
   GIT_REPOSITORY https://github.com/aradi/fypp
   GIT_TAG 3.1
)
FetchContent_MakeAvailable(fypp)
