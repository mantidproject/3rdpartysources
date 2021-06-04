set(cylabs_module_name "cylab")
set(cylabs_main ${CMAKE_CURRENT_SOURCE_DIR}/cylabs_main.f90)
set(cylabs_only_arg "cylabs")
list(APPEND cylabs_other_sources 
${CMAKE_CURRENT_SOURCE_DIR}/cylabs_subs.f90
)
f2py_fortran_library(${cylabs_module_name} ${cylabs_main} ${cylabs_only_arg} ${cylabs_other_sources})

set(fltabs_module_name "fltab")
set(fltabs_main ${CMAKE_CURRENT_SOURCE_DIR}/fltabs_main.f90)
set(fltabs_only_arg "fltabs")
list(APPEND fltabs_other_sources 
${CMAKE_CURRENT_SOURCE_DIR}/fltabs_subs.f90
)
f2py_fortran_library(${fltabs_module_name} ${fltabs_main} ${fltabs_only_arg} ${fltabs_other_sources})

set(muscat_module_name "muscat")
set(muscat_main ${CMAKE_CURRENT_SOURCE_DIR}/muscat_data_main.f90)
set(muscat_only_arg "muscat_data")
list(APPEND muscat_other_sources 
${CMAKE_CURRENT_SOURCE_DIR}/muscat_data.f90
${CMAKE_CURRENT_SOURCE_DIR}/muscat_geom.f90
)
f2py_fortran_library(${muscat_module_name} ${muscat_main} ${muscat_only_arg} ${muscat_other_sources})



