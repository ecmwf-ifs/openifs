# (C) Copyright 1989- ECMWF.
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# 
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction

if( NOT TARGET ifsdummy )

    ecbuild_info("[dummy]")

    ecbuild_list_add_pattern(LIST ifsdummy_src
        GLOB  dummy/ifsaux/*
              dummy/wam/*
              dummy/arpifs/*
    )

    if( HAVE_FORECAST_ONLY )
      ecbuild_list_add_pattern(LIST ifsdummy_src GLOB dummy/forecast/*)
    endif()

    ecbuild_add_library(TARGET ifsdummy
        SOURCES ${ifsdummy_src}
    )

    include(target_ignore_missing_symbols)
    target_ignore_missing_symbols( TARGET ifsdummy SYMBOLS
      abor1_
    )
  
endif()
