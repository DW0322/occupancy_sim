cmake_minimum_required(VERSION 3.16...4.0)
project(occupancy_sim VERSION 5.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

option(OCCUPANCY_BUILD_TESTS "Build the test suite" ON)
option(OCCUPANCY_FETCH_MORPHIO "Fetch and build MorphIO if not installed" ON)

# compile_commands.json for clangd. Set here rather than left to the caller
# because MorphIO turns it on itself - without this the file existed but only
# described the fetched dependency, so editors had no flags for our own code.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# HDF5 headers are a hard requirement - MorphIO reads .h5 through HighFive and
# won't configure without them. Checked up front because otherwise it surfaces
# as an unrelated missing-target error.
find_package(HDF5 QUIET COMPONENTS C)
if(NOT HDF5_FOUND)
  message(FATAL_ERROR
    "HDF5 development headers not found. Install the DEVELOPMENT package:\n"
    "  Debian/Ubuntu  apt install libhdf5-dev\n"
    "  Fedora/RHEL    dnf install hdf5-devel\n"
    "  Arch           pacman -S hdf5\n"
    "  macOS          brew install hdf5")
endif()

include(cmake/GetMorphIO.cmake)

add_library(occupancy_core
  src/segment.cpp
  src/schedule.cpp
  src/distance_transform.cpp
  src/migration.cpp
  src/simulation.cpp
  src/compartments.cpp
  src/stiffness.cpp
  src/scaffold_geometry.cpp
  src/clustering.cpp)
add_library(occupancy::core ALIAS occupancy_core)

target_include_directories(occupancy_core PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)
target_link_libraries(occupancy_core PUBLIC occupancy::morphio)

# Warnings for everything we compile. This used to cover the library only,
# which left the 1000-line CLI unchecked - it happens to be clean, but nothing
# was keeping it that way. Not -Werror: a researcher on a newer compiler should
# get a warning, not a build that refuses to produce the engine.
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  set(OCCUPANCY_WARNINGS -Wall -Wextra)
else()
  set(OCCUPANCY_WARNINGS)
endif()
target_compile_options(occupancy_core PRIVATE ${OCCUPANCY_WARNINGS})

add_executable(occupancy_engine apps/occupancy_engine.cpp)
target_link_libraries(occupancy_engine PRIVATE occupancy::core)
target_compile_options(occupancy_engine PRIVATE ${OCCUPANCY_WARNINGS})

# morphometry over the engine's own loader; where the territory radius R comes
# from. See apps/morph_stats.cpp.
add_executable(morph_stats apps/morph_stats.cpp)
target_link_libraries(morph_stats PRIVATE occupancy::core)
target_compile_options(morph_stats PRIVATE ${OCCUPANCY_WARNINGS})

if(OCCUPANCY_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

include(GNUInstallDirs)
install(TARGETS occupancy_engine RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(TARGETS occupancy_core
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY include/occupancy DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  FILES_MATCHING PATTERN "*.hpp")

message(STATUS "occupancy_sim ${PROJECT_VERSION}  build=${CMAKE_BUILD_TYPE}  tests=${OCCUPANCY_BUILD_TESTS}")
