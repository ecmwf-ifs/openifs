
# Activate warnings, ecbuild macros check the compiler recognises the options
if(HAVE_WARNINGS)
  ecbuild_add_c_flags("-Wall"                                  NO_FAIL)
  ecbuild_add_c_flags("-Wextra"                                NO_FAIL)
  ecbuild_add_c_flags("-Wno-unused-parameter"                  NO_FAIL)
  ecbuild_add_c_flags("-Wno-unused-variable"                   NO_FAIL)
  ecbuild_add_c_flags("-Wno-gnu-zero-variadic-macro-arguments" NO_FAIL)
endif()

# Always disable some warnings
ecbuild_add_c_flags("-Wno-deprecated-declarations" NO_FAIL)
if( CMAKE_C_COMPILER_ID MATCHES Intel )
  ecbuild_add_c_flags("-diag-disable=279")   # controlling expression is constant
  ecbuild_add_c_flags("-diag-disable=11076") # inline limits
endif()
