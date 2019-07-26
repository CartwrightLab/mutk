# Copyright (c) 2019 Reed A. Cartwright <reed@cartwright.ht>
# 
# This file is part of the Ultimate Source Code Project.
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

set(mutk_sources version.cc)

add_library(libmutk STATIC ${mutk_sources})
add_dependencies(libmutk configure-version.h)
target_compile_features(libmutk PUBLIC cxx_std_17)
target_link_libraries(libmutk PUBLIC doctest::doctest)
target_compile_definitions(libmutk PRIVATE DOCTEST_CONFIG_DISABLE)

target_include_directories(libmutk PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../include")
target_include_directories(libmutk PUBLIC "${CMAKE_CURRENT_BINARY_DIR}/..")
set_target_properties(libmutk PROPERTIES OUTPUT_NAME mutk)

add_library(libmutk-doctest OBJECT ${mutk_sources})
get_property(doctest_include_dir TARGET doctest::doctest
    PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
target_include_directories(libmutk-doctest PRIVATE ${doctest_include_dir})
get_property(libmutk_include_dir TARGET libmutk
    PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
target_include_directories(libmutk-doctest PRIVATE ${libmutk_include_dir})
