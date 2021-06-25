#.rst:
#
# Find the gfortran executable and its version

find_program(GFORTRAN_EXECUTABLE NAMES gfortran)

if(GFORTRAN_EXECUTABLE)
execute_process(COMMAND "${GFORTRAN_EXECUTABLE}" -dumpversion
OUTPUT_VARIABLE GFORTRAN_VERSION_STRING
OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GFortran
  REQUIRED_VARS GFORTRAN_EXECUTABLE
  VERSION_VAR GFORTRAN_VERSION_STRING
  )

mark_as_advanced(GFORTRAN_EXECUTABLE GFORTRAN_VERSION_STRING)