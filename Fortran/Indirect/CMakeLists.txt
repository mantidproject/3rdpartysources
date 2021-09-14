

set("3rdparty_libs_output_dir" ${PROJECT_BINARY_DIR}/bin)

add_custom_target(build-time-make-directory ALL
    COMMAND ${CMAKE_COMMAND} -E make_directory ${3rdparty_libs_output_dir})

function(f2py_fortran_library module_name module_main only_arg) 
set(generated_module_file ${module_name}${PYTHON_EXTENSION_MODULE_SUFFIX})
add_custom_target(${module_name} ALL
  DEPENDS ${3rdparty_libs_output_dir}/${generated_module_file}
)
# windows commands
if(WIN32)
add_custom_command(
  OUTPUT ${3rdparty_libs_output_dir}/${generated_module_file}
  COMMAND ${Python_EXECUTABLE} -m numpy.f2py -c --fcompiler=gnu95 --compiler=mingw32 --opt=${GFORTRAN_FLAGS} -m ${module_name} ${module_main} only: ${only_arg} : ${ARGN}
  WORKING_DIRECTORY ${3rdparty_libs_output_dir}
)
# unix commands
else()
add_custom_command(
  OUTPUT ${3rdparty_libs_output_dir}/${generated_module_file}
  COMMAND ${Python_EXECUTABLE} -m numpy.f2py -c --fcompiler=gnu95 --compiler=unix --opt=${GFORTRAN_FLAGS} -m ${module_name} ${module_main} only: ${only_arg} : ${ARGN}
  WORKING_DIRECTORY ${3rdparty_libs_output_dir}
)
endif()
endfunction()

add_subdirectory(mantidindirect/Bayes)
add_subdirectory(mantidindirect/AbsCorrection)
