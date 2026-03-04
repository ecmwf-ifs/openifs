# (C) Copyright 2019 ECMWF.
#
# This file is covered by the LICENSING file in the root of this project.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

function( run_command COMMAND )
  set( COMMAND ${ARGV} )
  string(REPLACE ";" " " PRINT_COMMAND "${COMMAND}" )
  message( "${PRINT_COMMAND}" )
  execute_process(
      COMMAND ${COMMAND}
      RESULT_VARIABLE res
      OUTPUT_VARIABLE stdout
      ERROR_VARIABLE  stderr
  )
  message("${stdout}")
  message("${stderr}")

  string(REGEX MATCH "${PASS_REGULAR_EXPRESSION}" stdout_match "${stdout}" )
  string(REGEX MATCH "${PASS_REGULAR_EXPRESSION}" stderr_match "${stderr}" )

  if( stdout_match OR stderr_match )
    message( "Test succeeded: Regex [${PASS_REGULAR_EXPRESSION}] was found in program output")
  else()
    message( FATAL_ERROR "Test failed: Could not find regex [${PASS_REGULAR_EXPRESSION}] in program output")
  endif()
endfunction()

message( "Running test ${EXECUTABLE} ... ")
run_command( ${LAUNCH} ${EXECUTABLE} )
message( "Running test ${EXECUTABLE} ... done")
