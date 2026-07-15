# (C) Copyright 1989- ECMWF.
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# 
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction

ecbuild_info("[wam]")

##################################################################################################################
### Add ecwam project depending on ifs-source configuration

set( ECWAM_PROJECT_NAME ecwam_${prec} )
ifs_propagate_flags( ${ECWAM_PROJECT_NAME} )

set( ECWAM_${PREC}_ENABLE_SINGLE_PRECISION OFF  )
if( HAVE_SINGLE_PRECISION )
  set( ECWAM_${PREC}_ENABLE_SINGLE_PRECISION ON  )
endif()
set( ECWAM_${PREC}_OCEANMODEL_LIBRARIES ${NEMOGCMCOUP_LIBRARIES} )
set( ECWAM_${PREC}_OCEANMODEL_HAVE_SINGLE_PRECISION ${NEMO_${PREC}_HAVE_SINGLE_PRECISION} )
set( ECWAM_${PREC}_ENABLE_UNWAM ON  )
set( ECWAM_${PREC}_ENABLE_ECFLOW ON )

add_subdirectory( contrib/ecwam )
ecbuild_find_package( ${ECWAM_PROJECT_NAME} REQUIRED )

##################################################################################################################
### Files that need to be part of arpifs.${PREC} target

list(APPEND wam_ifs_srcs  # Callback functions for handlers in coupled mode
  wam/Wam_oper/ifstowam.F90
  wam/Wam_oper/outwspec_io_serv.F90
  wam/Wam_oper/outint_io_serv.F90
)
set_property(SOURCE ${wam_ifs_srcs} PROPERTY COMPILE_OPTIONS ${autopromote_flags})

##################################################################################################################
### Generate interfaces for wam_ifs_srcs

# First create symlinks of above files to a source directory in build-dir.
# These files only will then be considered for generating interfaces.

set( WAM_SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/wam/wam-src )
file( MAKE_DIRECTORY ${WAM_SOURCE_DIR})
set( CREATE_SYMLINKS )
foreach( f ${wam_ifs_srcs} )
  get_filename_component(basename "${f}" NAME)
  set( symlink ${WAM_SOURCE_DIR}/${basename} )
  set( CREATE_SYMLINK COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/${f} ${symlink} )
  list( APPEND CREATE_SYMLINKS ${CREATE_SYMLINK} )
endforeach()
execute_process( ${CREATE_SYMLINKS} )


if( NOT TARGET wam_intfb )
  ecbuild_generate_fortran_interfaces(
    TARGET wam_intfb
    DIRECTORIES wam-src
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/wam
    DESTINATION wam
    INCLUDE_DIRS wam_intfb_includes
    PARALLEL ${FCM_PARALLEL}
  )
endif()

set( ECWAM_LIBRARIES ecwam_${prec} wam_intfb )

##################################################################################################################
### Tools that don't need data assimilation, and not yet part of ecwam

if( ${ECWAM_PROJECT_NAME}_HAVE_UNWAM )
  list( APPEND WAM_DEFINITIONS WAM_HAVE_UNWAM )
  ecbuild_info("ecwam was compiled with UNWAM support")
else()
  ecbuild_info("ecwam was not compiled with UNWAM support")
endif()

foreach(program IN ITEMS

    bouint
    preproc
    preset
    create_wam_bathymetry
    create_wam_bathymetry_ETOPO1 )

    if( NOT TARGET ${program}-create_symlink )
      add_custom_target(${program}-create_symlink
        ALL
        COMMAND ${CMAKE_COMMAND} -E
        create_symlink $<TARGET_FILE_NAME:${ECWAM_PROJECT_NAME}-${program}> ${CMAKE_BINARY_DIR}/bin/${program}
      )
      install(FILES ${CMAKE_BINARY_DIR}/bin/${program} TYPE BIN)
    endif()

endforeach()

foreach(name IN ITEMS

    intwaminput
    write_mpdecomp )

    if( NOT TARGET ${name} )
        ecbuild_add_executable(TARGET ${name}
            SOURCES wam/Wam_oper/${name}.F90 wam/Alt/getclo.F
            LIBS
                ecwam_${prec}
                OpenMP::OpenMP_Fortran
            LINKER_LANGUAGE Fortran)
        set_property(SOURCE wam/Wam_oper/${name}.F90 PROPERTY COMPILE_OPTIONS ${autopromote_flags})
    endif()

endforeach()

##################################################################################################################

