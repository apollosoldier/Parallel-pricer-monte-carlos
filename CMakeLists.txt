cmake_minimum_required(VERSION 3.0)
project(mc-pricer CXX)

# Prefix
if (NOT PREFIX)
    set(PREFIX ${CMAKE_CURRENT_BINARY_DIR})
endif (NOT PREFIX)
set(CMAKE_INSTALL_PREFIX ${PREFIX} CACHE STRING "Installation prefix." FORCE)

if (CXX)
    set(CMAKE_CXX_COMPILER ${CXX})
endif (CXX)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

# Release or Debug
if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "Setting build type to 'Debug' as none was specified.")
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build." FORCE)
endif ()

# Detect OPM
find_package(OpenMP QUIET)
if (OPENMP_FOUND)
    set (CMAKE_CXX_FLAGS ${CMAKE_CXX_CFLAGS} ${OpenMP_CXX_FLAGS})
    message(STATUS "OpenMP found.")
else (OPENMP_FOUND)
    message(STATUS "OpenMP not found.")
endif (OPENMP_FOUND)

# Detect MPI
add_definitions(-DOMPI_SKIP_MPICXX)
set(MPI_LIBS)
find_package(MPI REQUIRED)
if (MPI_C_FOUND)
    include_directories(${MPI_C_INCLUDE_PATH})
    set (MPI_LIBS ${MPI_LIBS} ${MPI_C_LIBRARIES})
else (MPI_C_FOUND)
    message(STATUS "MPI not found.")
endif (MPI_C_FOUND)

# Detect PNL
set(CMAKE_PREFIX_PATH /matieres/5MMPMP6/pnl)

find_package(Pnl REQUIRED)
set(LIBS ${LIBS} ${PNL_LIBRARIES})
include_directories(${PNL_INCLUDE_DIRS})
if(${CMAKE_BUILD_TYPE} STREQUAL "Release")
    add_definitions(-DPNL_RANGE_CHECK_OFF)
endif()

# Add math lib
set(LIBS ${LIBS} m)

# Add jlparser
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
add_subdirectory(src/3rdparty/jlparser)

#Simulation
add_executable(pricer src/pricer.cpp src/BlackScholesModel.cpp src/VanillaCall.cpp src/VanillaCall.hpp src/BasketOption.cpp src/BasketOption.hpp src/AsianOption.cpp src/AsianOption.hpp src/PerfOption.cpp src/PerfOption.hpp src/MonteCarlo.cpp src/MonteCarlo.hpp )
target_link_libraries(pricer jl_parser ${LIBS})
# Print compiler flags
get_directory_property(DirDefs COMPILE_DEFINITIONS)
message(STATUS "COMPILE_DEFINITIONS = ${DirDefs}")
if (${CMAKE_BUILD_TYPE} STREQUAL "Release")
    message(STATUS "C++ Flags: ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE} ${CMAKE_LINKER_EXE_FLAGS}" )
else()
    message(STATUS "C++ Flags: ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_LINKER_EXE_FLAGS}" )
endif()
