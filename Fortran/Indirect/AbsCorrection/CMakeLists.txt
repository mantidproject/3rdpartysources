set(muscat_module_name "muscat")
set(muscat_main ${CMAKE_CURRENT_SOURCE_DIR}/muscat_data_main.f90)
set(muscat_only_arg "muscat_data")
list(APPEND muscat_other_sources 
${CMAKE_CURRENT_SOURCE_DIR}/muscat_data.f90
${CMAKE_CURRENT_SOURCE_DIR}/muscat_geom.f90
)
f2py_fortran_library(${muscat_module_name} ${muscat_main} ${muscat_only_arg} ${muscat_other_sources})



