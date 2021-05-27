file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/bin)


function(f2py_fortran_library module_name module_main only_arg) 
set(generated_module_file ${CMAKE_CURRENT_BINARY_DIR}/${module_name}${PYTHON_EXTENSION_MODULE_SUFFIX})
add_custom_target(${module_name} ALL
  DEPENDS ${generated_module_file}
)
if(WIN32)
add_custom_command(
  OUTPUT ${generated_module_file}
  COMMAND ${F2PY_EXECUTABLE} -c --fcompiler=gnu95 --compiler=mingw32 -m ${module_name} ${module_main} only: ${only_arg} : ${ARGN}
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)
else()
add_custom_command(
  OUTPUT ${generated_module_file}
  COMMAND ${F2PY_EXECUTABLE} -c --fcompiler=gnu95 --compiler=unix -m ${module_name} ${module_main} only: ${only_arg} : ${ARGN}
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)
endif()
endfunction()

add_subdirectory(Bayes)
add_subdirectory(AbsCorrection)
