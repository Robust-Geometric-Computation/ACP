cmake_minimum_required(VERSION 3.11)

# Project options
project(acp VERSION 1.0
            DESCRIPTION "Adaptive Controlled Perturbation"
            LANGUAGES CXX C)

set(PROJECT_VERSION 1.0)
set(PROJECT_VERSION_MAJOR 1)
set(PROJECT_VERSION_MINOR 0)

if(UNIX AND NOT APPLE)
  set(LINUX TRUE)
endif()

# Only do these if this is the main project, and not if it is included through add_subdirectory
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)

    # Set default build type to debug 
    set(default_build_type "Debug")
    if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)

      message(STATUS "Setting build type to '${default_build_type}' as none was specified.")

      set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
          STRING "Choose the type of build." FORCE)

      # Set the possible values of build type for cmake-gui
      set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
    endif()
    
    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)
    
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
    
    # Die on first error
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wfatal-errors")
    
    set(CMAKE_MODULE_PATH /usr/local/lib/cmake /usr/local/lib/x86_64-linux-gnu/cmake ${CMAKE_SOURCE_DIR}/cmake/modules)
    set(CMAKE_PREFIX_PATH /usr/local/lib/cmake/glfw)

    # Optionally set things like CMAKE_CXX_STANDARD, CMAKE_POSITION_INDEPENDENT_CODE here
    # add_compile_options(-Wall -Wextra -pedantic)

    # Let's nicely support folders in IDE's
    set_property(GLOBAL PROPERTY USE_FOLDERS ON)

    # Testing only available if this is the main app
    # Note this needs to be done in the main CMakeLists
    # since it calls enable_testing, which must be in the
    # main CMakeLists.
    option(PACKAGE_TESTS "Build the tests" OFF)
    if(PACKAGE_TESTS)
        enable_testing()
        include(GoogleTest)
        add_subdirectory(tests)
    endif()

    # Docs only available if this is the main app
    #find_package(Doxygen)
    #if(Doxygen_FOUND)
	  #    add_subdirectory(docs)
    #else()
	  #    message(STATUS "Doxygen not found, not building docs")
    #endif()
endif()

# Find all dependencies

find_package (PkgConfig REQUIRED)
find_package (OpenGL REQUIRED)
find_package (glfw3 REQUIRED)
find_package (GLM REQUIRED)
find_package (GLEW REQUIRED STATIC)
find_package (GLUT REQUIRED)

find_package (gflags REQUIRED)

# See if VTK is using OpenGL or OpenGL2
find_package(VTK CONFIG)
if (NOT VTK_RENDERING_BACKEND)
  set(VTK_RENDERING_BACKEND OpenGL)
endif()
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})


if (APPLE)
  include_directories(/usr/local/include)
  find_library(COCOA_LIBRARY Cocoa REQUIRED)
  find_library(IOKIT_LIBRARY IOKit REQUIRED)
  find_library(COREVID_LIBRARY CoreVideo REQUIRED)
  add_definitions("-DOSX")
endif (APPLE)
  
#add_definitions("-DVERBOSE")
add_definitions("-DGLM_ENABLE_EXPERIMENTAL")
add_definitions("-DNO_HOMOTOPY")
add_definitions("-DNO_MODE")

find_package (MPFR REQUIRED)
find_package (GMP REQUIRED)

find_package(LAPACK REQUIRED)
find_package(Freetype REQUIRED)
find_package(PNG REQUIRED)

find_package (Threads REQUIRED)

set(GLAD "${CMAKE_SOURCE_DIR}/extern/glad/glad.h"
         "${CMAKE_SOURCE_DIR}/extern/glad/glad.c")

# Submodules 
# include(FetchContent)

# Subdirectories
include_directories( ${CMAKE_SOURCE_DIR}/include )
include_directories( ${CMAKE_SOURCE_DIR}/extern/ )
include_directories( ${FREETYPE_INCLUDE_DIRS} )
include_directories( ${VTK_INCLUDE_DIRS} )

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

add_subdirectory(src)
add_subdirectory(apps)

install(DIRECTORY extras DESTINATION ${INSTALL_PATH}/extras)
