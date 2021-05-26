function(bayes_fortran_library module_name module_main only_arg) 
set(generated_module_file ${CMAKE_CURRENT_BINARY_DIR}/${module_name}${PYTHON_EXTENSION_MODULE_SUFFIX})
add_custom_target(${module_name} ALL
  DEPENDS ${generated_module_file}
)
if(WIN32)
add_custom_command(
  OUTPUT ${generated_module_file}
  COMMAND ${F2PY_EXECUTABLE} -c --fcompiler=gnu95 --compiler=mingw32 -m ${module_name} ${module_main} only: ${only_arg} : ${ARGN}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
else()
add_custom_command(
  OUTPUT ${generated_module_file}
  COMMAND ${F2PY_EXECUTABLE} -c --fcompiler=gnu95 --compiler=unix -m ${module_name} ${module_main} only: ${only_arg} : ${ARGN}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
endif()
endfunction()

set(Qldata_module_name "QLdata")
set(QLdata_main ${CMAKE_CURRENT_SOURCE_DIR}/Bayes/QLdata_main.f90)
set(QLdata_only_arg "qldata")
list(APPEND QLdata_other_sources 
${CMAKE_CURRENT_SOURCE_DIR}/Bayes/QLdata_subs.f90
${CMAKE_CURRENT_SOURCE_DIR}/Bayes/Util.f90
${CMAKE_CURRENT_SOURCE_DIR}/Bayes/Bayes.f90 
${CMAKE_CURRENT_SOURCE_DIR}/Bayes/Four.f90)
bayes_fortran_library(${Qldata_module_name} ${QLdata_main} ${QLdata_only_arg} ${QLdata_other_sources})


set(QLres_module_name "QLres")
set(QLres_main ${CMAKE_CURRENT_SOURCE_DIR}/Bayes/QLres_main.f90)
set(QLres_only_arg "qlres")
list(APPEND QLres_other_sources 
${CMAKE_CURRENT_SOURCE_DIR}/Bayes/QLres_subs.f90
${CMAKE_CURRENT_SOURCE_DIR}/Bayes/BlrRes.f90
${CMAKE_CURRENT_SOURCE_DIR}/Bayes/Bayes.f90 
${CMAKE_CURRENT_SOURCE_DIR}/Bayes/Four.f90 
${CMAKE_CURRENT_SOURCE_DIR}/Bayes/Util.f90
)
bayes_fortran_library(${QLres_module_name} ${QLres_main} ${QLres_only_arg} ${QLres_other_sources})

set(QLse_module_name "QLse")
set(QLse_main ${CMAKE_CURRENT_SOURCE_DIR}/Bayes/QLse_main.f90)
set(QLse_only_arg "qlstexp")
list(APPEND QLse_other_sources
${CMAKE_CURRENT_SOURCE_DIR}/Bayes/QLse_subs.f90
${CMAKE_CURRENT_SOURCE_DIR}/Bayes/BlrRes.f90
${CMAKE_CURRENT_SOURCE_DIR}/Bayes/Bayes.f90
${CMAKE_CURRENT_SOURCE_DIR}/Bayes/Four.f90
${CMAKE_CURRENT_SOURCE_DIR}/Bayes/Util.f90
${CMAKE_CURRENT_SOURCE_DIR}/Bayes/Simopt.f90
)
bayes_fortran_library(${QLse_module_name} ${QLse_main} ${QLse_only_arg} ${QLse_other_sources})

set(Quest_module_name "Quest")
set(Quest_main ${CMAKE_CURRENT_SOURCE_DIR}/Bayes/Quest_main.f90)
set(Quest_only_arg "quest")
list(APPEND Quest_other_sources 
${CMAKE_CURRENT_SOURCE_DIR}/Bayes/Quest_subs.f90
${CMAKE_CURRENT_SOURCE_DIR}/Bayes/BlrRes.f90
${CMAKE_CURRENT_SOURCE_DIR}/Bayes/Bayes.f90
${CMAKE_CURRENT_SOURCE_DIR}/Bayes/Four.f90 
${CMAKE_CURRENT_SOURCE_DIR}/Bayes/Util.f90
${CMAKE_CURRENT_SOURCE_DIR}/Bayes/Simopt.f90
)
bayes_fortran_library(${Quest_module_name} ${Quest_main} ${Quest_only_arg} ${Quest_other_sources})


set(ResNorm_module_name "ResNorm")
set(ResNorm_main ${CMAKE_CURRENT_SOURCE_DIR}/Bayes/ResNorm_main.f90)
set(ResNorm_only_arg "resnorm")
list(APPEND ResNorm_other_sources 
${CMAKE_CURRENT_SOURCE_DIR}/Bayes/ResNorm_subs.f90
${CMAKE_CURRENT_SOURCE_DIR}/Bayes/BlrRes.f90
${CMAKE_CURRENT_SOURCE_DIR}/Bayes/Bayes.f90
${CMAKE_CURRENT_SOURCE_DIR}/Bayes/Four.f90 
${CMAKE_CURRENT_SOURCE_DIR}/Bayes/Util.f90
)
bayes_fortran_library(${ResNorm_module_name} ${ResNorm_main} ${ResNorm_only_arg} ${ResNorm_other_sources})