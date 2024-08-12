cmake_minimum_required(VERSION 3.15...3.27)
project(stark)
set(CMAKE_CXX_STANDARD 20)

# Compile options for address sanitizer and Libcxx debug
#add_compile_options(-fsanitize=address -fno-omit-frame-pointer -g -D_GLIBCXX_DEBUG)
#add_link_options(-fsanitize=address)

add_subdirectory(stark)
add_subdirectory(examples)
add_subdirectory(tests)

option(STARK_BUILD_PYTHON_BINDINGS "Build Python Bindings" ON)
if(STARK_BUILD_PYTHON_BINDINGS)
	add_subdirectory(pystark/cpp)
endif()
