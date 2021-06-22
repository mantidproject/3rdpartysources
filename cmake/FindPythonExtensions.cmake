#.rst:
#
# This module defines CMake variables to build Python extension modules. 
# Adapted from https://github.com/scikit-build/scikit-build/blob/master/skbuild/resources/cmake/FindPythonExtensions.cmake
#

#=============================================================================
# Copyright 2011 Kitware, Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#=============================================================================
set(_command "
import distutils.sysconfig
import itertools
import os
import os.path
import site
import sys

result = None
rel_result = None
candidate_lists = []

try:
    candidate_lists.append((distutils.sysconfig.get_python_lib(),))
except AttributeError: pass

try:
    candidate_lists.append(site.getsitepackages())
except AttributeError: pass

try:
    candidate_lists.append((site.getusersitepackages(),))
except AttributeError: pass

candidates = itertools.chain.from_iterable(candidate_lists)

for candidate in candidates:
    rel_candidate = os.path.relpath(
      candidate, sys.prefix)
    if not rel_candidate.startswith(\"..\"):
        result = candidate
        rel_result = rel_candidate
        break

ext_suffix_var = 'SO'
if sys.version_info[:2] >= (3, 5):
    ext_suffix_var = 'EXT_SUFFIX'

sys.stdout.write(\";\".join((
    os.sep,
    os.pathsep,
    sys.prefix,
    result,
    rel_result,
    distutils.sysconfig.get_config_var(ext_suffix_var)
)))
")

if(Python_Interpreter_FOUND)
execute_process(COMMAND "${Python_EXECUTABLE}" -c "${_command}"
                OUTPUT_VARIABLE _list
                RESULT_VARIABLE _result)

list(GET _list 0 _item)
set(PYTHON_SEPARATOR "${_item}")
mark_as_advanced(PYTHON_SEPARATOR)

list(GET _list 1 _item)
set(PYTHON_PATH_SEPARATOR "${_item}")
mark_as_advanced(PYTHON_PATH_SEPARATOR)

list(GET _list 2 _item)
set(PYTHON_PREFIX "${_item}")
mark_as_advanced(PYTHON_PREFIX)

list(GET _list 3 _item)
set(PYTHON_SITE_PACKAGES_DIR "${_item}")
mark_as_advanced(PYTHON_SITE_PACKAGES_DIR)

list(GET _list 4 _item)
set(PYTHON_RELATIVE_SITE_PACKAGES_DIR "${_item}")
mark_as_advanced(PYTHON_RELATIVE_SITE_PACKAGES_DIR)

if(NOT DEFINED PYTHON_EXTENSION_MODULE_SUFFIX)
  list(GET _list 5 _item)
  set(PYTHON_EXTENSION_MODULE_SUFFIX "${_item}")
endif()

else()
message(STATUS "To define Python exetensions Python interpretator is required to be found.")
endif()