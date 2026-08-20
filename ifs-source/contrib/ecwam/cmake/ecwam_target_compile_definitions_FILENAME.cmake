# (C) Copyright 2022- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

function( ecwam_target_compile_definitions_FILENAME target )
    get_target_property( sources "${target}" SOURCES)
    foreach( src ${sources} )
        get_filename_component( filename ${src} NAME )
        set_property(
            SOURCE ${src}
            APPEND
            PROPERTY COMPILE_DEFINITIONS __FILENAME__="${filename}" )
    endforeach()
endfunction()

