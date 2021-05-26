cmake_minimum_required(VERSION 3.16)

project(MantidFortranRoutines)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

find_package (Python COMPONENTS Interpreter Development NumPy)

# Find numpy headers
exec_program(${Python_EXECUTABLE}
             ARGS "-c \"import numpy; print(numpy.get_include())\""
             OUTPUT_VARIABLE NUMPY_INCLUDE_DIR
             RETURN_VALUE NUMPY_NOT_FOUND
            )
if(NUMPY_NOT_FOUND)
    message(FATAL_ERROR "NumPy headers not found")
endif()

find_package(F2PY REQUIRED)
# include python and numpy dirs
include_directories(${Python_INCLUDE_DIRS})
include_directories(${NUMPY_INCLUDE_DIR})

add_subdirectory(Fortran)

