cmake_minimum_required(VERSION 3.7.2)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMake/Modules/")
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/CMake)
set(CMAKE_BUILD_TYPE Release)
project(el3_s3d CXX)

# My optimization flags
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_OPT_ARCH_NATIVE_SUPPORTED)
if (COMPILER_OPT_ARCH_NATIVE_SUPPORTED)
    set(CMAKE_CXX_FLAGS_RELEASE "-march=native -DNDEBUG")
else(COMPILER_OPT_ARCH_NATIVE_SUPPORTED)
	set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG")
endif(COMPILER_OPT_ARCH_NATIVE_SUPPORTED)
CHECK_CXX_COMPILER_FLAG("-Ofast" COMPILER_OPT_OFAST_SUPPORTED)
if (COMPILER_OPT_OFAST_SUPPORTED)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Ofast")
else(COMPILER_OPT_ARCH_NATIVE_SUPPORTED)
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
endif(COMPILER_OPT_OFAST_SUPPORTED)

# VTK
find_package(VTK)
if (VTK_FOUND)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_VTK")
   include(${VTK_USE_FILE})
   link_libraries(${VTK_LIBRARIES})
endif(VTK_FOUND)

# Eigen3
#set (EIGEN_USE_MKL_ALL) 
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# BLAS
find_package(BLAS REQUIRED)
include_directories(${BLAS_INCLUDE_DIR})
link_libraries(${BLAS_LIBRARIES})

# LAPACKE
find_package(LAPACK REQUIRED)
include_directories(${LAPACK_INCLUDE_DIR})
link_libraries(${LAPACK_LIBRARIES})

# Boost
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
link_libraries(${Boost_LIBRARIES})

# OpenMP
find_package(OpenMP)
if (OpenMP_FOUND)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
link_libraries(${OpenMP_CXX_LIBRARIES}) 
endif(OpenMP_FOUND)

add_executable(el3_s3d main.cpp mesh.hpp pde.cpp export.cpp tet.hpp utils.hpp face.hpp mesh.cpp names.h pde.hpp tmv.hpp vertex.hpp)
