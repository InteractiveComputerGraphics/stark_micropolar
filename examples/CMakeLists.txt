set(EIGEN3_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../stark/extern/Eigen)
if( NOT EIGEN3_INCLUDE_DIR )
    message( FATAL_ERROR "Please point the environment variable EIGEN3_INCLUDE_DIR to the include directory of your Eigen3 installation.")
endif()
set (CMAKE_CXX_FLAGS "-DEIGEN_MPL2_ONLY")

# Enable SIMD and suppress compiler warnings
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
	add_compile_options("-march=native")
    add_compile_options("-Wno-ignored-attributes")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
	add_compile_options("/arch:AVX")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
	add_compile_options("-march=native")
endif()

set(SOURCE_FILES
	main.cpp
	)

# Create executable    
add_executable(examples ${SOURCE_FILES})

# Include/link dependencies

## MKL
if (${STARK_ENABLE_MKL})
	target_compile_options(examples PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
	target_include_directories(examples PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
	target_link_libraries(examples $<LINK_ONLY:MKL::MKL>)
endif()

## Stark
target_link_libraries(examples stark)
target_include_directories(examples PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../stark/include)

## Eigen
target_include_directories(examples PUBLIC ${EIGEN3_INCLUDE_DIR})
