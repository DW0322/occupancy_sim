# GetMorphIO.cmake - provides the target `occupancy::morphio`.
#
# OCCUPANCY_MORPHIO_SOURCE picks where MorphIO comes from:
#   auto    (default) use a valid system install if there is one, else fetch
#   fetch             always fetch and build at the pinned tag
#   system            require a system install, fail if absent
#
# This file is defensive for four reasons, all of which bit us:
#
# 1. find_package(MorphIO) isn't safe. A stale install leaves a Config.cmake
#    pointing at a library that's gone, and that config FATAL_ERRORs from inside
#    find_package - QUIET doesn't suppress it and CMake can't catch it, so the
#    configure just dies. So the system path uses find_path/find_library only
#    and never runs MorphIO's config script.
# 2. There is no `morphio` target. MorphIO defines morphio_static and
#    morphio_shared; the namespaced MorphIO::morphio only exists through its
#    install/export config, not under FetchContent.
# 3. GIT_TAG master isn't reproducible, so it's pinned below. Bump deliberately.
# 4. Newer toolchains break the vendored sources: CMake 4 rejects their
#    cmake_minimum_required(VERSION 3.2) (hence CMAKE_POLICY_VERSION_MINIMUM),
#    GCC 15/16 dropped the transitive <cstdint> MorphIO needs for uint32_t
#    (hence the forced include), and gcc 13+ throws a spurious warning in the
#    bundled lexertl14 that MorphIO's -Werror makes fatal (hence WERROR=OFF).
#
# Static link when fetching, on purpose: the shared library needs
# LD_LIBRARY_PATH at runtime and forgetting it fails in confusing ways.
#
# No network? FetchContent's own escape hatch takes a local checkout and skips
# the download:
#   cmake --preset default -DFETCHCONTENT_SOURCE_DIR_MORPHIO=/path/to/MorphIO
# Worth knowing because upstream is archived read-only, so this pin can never
# move forward and a vanished repository has no upstream fix. See
# docs/RUNNING.md section 1.

if(TARGET occupancy::morphio)
  return()
endif()

# v3.4.0 is the last release, not just our current pin - upstream was archived
# read-only on 2025-02-26, so there's no newer tag to bump to.
set(OCCUPANCY_MORPHIO_TAG "v3.4.0" CACHE STRING
    "MorphIO git tag to fetch. Pinned for reproducibility; also the last tag the archived upstream repo ever cut.")
set(OCCUPANCY_MORPHIO_SOURCE "auto" CACHE STRING
    "Where MorphIO comes from: auto, fetch or system")
set_property(CACHE OCCUPANCY_MORPHIO_SOURCE PROPERTY STRINGS auto fetch system)

# apply the toolchain workarounds to a MorphIO target we don't own
function(_occupancy_patch_morphio_target tgt)
  if(NOT TARGET ${tgt})
    return()
  endif()
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    # force in the <cstdint> GCC 15/16 no longer includes transitively, rather
    # than patching vendored sources
    target_compile_options(${tgt} PRIVATE -include cstdint)
  elseif(MSVC)
    target_compile_options(${tgt} PRIVATE /FI cstdint)
  endif()
endfunction()

# ---------------------------------------------------------------------------
# system: headers + library only, no script execution
set(_occ_have_system FALSE)
if(NOT OCCUPANCY_MORPHIO_SOURCE STREQUAL "fetch")
  find_path(MORPHIO_INCLUDE_DIR morphio/morphology.h
    HINTS ${MORPHIO_ROOT} /usr/local/include /usr/include)
  find_library(MORPHIO_LIBRARY NAMES morphio
    HINTS ${MORPHIO_ROOT} /usr/local/lib /usr/lib)

  if(MORPHIO_INCLUDE_DIR AND MORPHIO_LIBRARY AND EXISTS "${MORPHIO_LIBRARY}")
    # MorphIO's headers pull in HighFive, which pulls in HDF5
    find_path(MORPHIO_HIGHFIVE_DIR highfive/H5File.hpp
      HINTS ${MORPHIO_INCLUDE_DIR} /usr/local/include /usr/include)
    if(MORPHIO_HIGHFIVE_DIR)
      set(_occ_have_system TRUE)
    else()
      message(STATUS
        "MorphIO: found ${MORPHIO_LIBRARY} but its bundled HighFive headers "
        "are missing - treating the install as unusable and fetching instead")
    endif()
  endif()
endif()

if(_occ_have_system)
  message(STATUS "MorphIO: using system install ${MORPHIO_LIBRARY}")
  message(STATUS "MorphIO: shared system library - you may need "
                 "LD_LIBRARY_PATH at runtime. Configure with "
                 "-DOCCUPANCY_MORPHIO_SOURCE=fetch for a self-contained build.")
  add_library(occupancy_morphio INTERFACE)
  target_include_directories(occupancy_morphio SYSTEM INTERFACE
    ${MORPHIO_INCLUDE_DIR} ${MORPHIO_HIGHFIVE_DIR} ${HDF5_INCLUDE_DIRS})
  target_link_libraries(occupancy_morphio INTERFACE ${MORPHIO_LIBRARY} ${HDF5_LIBRARIES})
  add_library(occupancy::morphio ALIAS occupancy_morphio)
  return()
endif()

if(OCCUPANCY_MORPHIO_SOURCE STREQUAL "system")
  message(FATAL_ERROR
    "MorphIO not found and OCCUPANCY_MORPHIO_SOURCE=system.\n"
    "Either install MorphIO, or configure with "
    "-DOCCUPANCY_MORPHIO_SOURCE=fetch to build it as part of this project.")
endif()

# ---------------------------------------------------------------------------
# fetch
message(STATUS "MorphIO: fetching ${OCCUPANCY_MORPHIO_TAG} (a few minutes the first time)")
include(FetchContent)

set(BUILD_BINDINGS OFF CACHE BOOL "" FORCE)   # no python bindings, no pybind11
set(MORPHIO_TESTS OFF CACHE BOOL "" FORCE)
set(MorphIO_WERROR OFF CACHE BOOL "" FORCE)
set(MORPHIO_ENABLE_COVERAGE OFF CACHE BOOL "" FORCE)
set(HIGHFIVE_UNIT_TESTS OFF CACHE BOOL "" FORCE)
set(HIGHFIVE_EXAMPLES OFF CACHE BOOL "" FORCE)
set(HIGHFIVE_BUILD_DOCS OFF CACHE BOOL "" FORCE)
# MorphIO and HighFive still declare cmake_minimum_required(VERSION 3.2), which
# CMake 4 refuses outright.
set(CMAKE_POLICY_VERSION_MINIMUM 3.5 CACHE STRING "" FORCE)

FetchContent_Declare(morphio
  GIT_REPOSITORY https://github.com/BlueBrain/MorphIO.git
  GIT_TAG ${OCCUPANCY_MORPHIO_TAG}
  GIT_SHALLOW TRUE
  GIT_SUBMODULES_RECURSE TRUE
  GIT_PROGRESS TRUE)
FetchContent_MakeAvailable(morphio)

# Submodules carry the vendored HighFive, gsl-lite and lexertl. Without them you
# get a cascade of missing-target errors, so name the real cause.
if(NOT EXISTS "${morphio_SOURCE_DIR}/3rdparty/HighFive/CMakeLists.txt")
  message(FATAL_ERROR
    "MorphIO was fetched without its submodules. Delete the build directory "
    "and reconfigure; if it persists, check that git can fetch submodules "
    "through your network.")
endif()

foreach(_t morphio_obj morphio_static morphio_shared)
  _occupancy_patch_morphio_target(${_t})
endforeach()

# prefer static: makes the engine self-contained, no LD_LIBRARY_PATH needed
if(TARGET morphio_static)
  add_library(occupancy::morphio ALIAS morphio_static)
  message(STATUS "MorphIO: linking morphio_static (self-contained, no LD_LIBRARY_PATH needed)")
elseif(TARGET morphio_shared)
  add_library(occupancy::morphio ALIAS morphio_shared)
  message(STATUS "MorphIO: linking morphio_shared")
else()
  message(FATAL_ERROR
    "MorphIO was fetched but defines neither morphio_static nor "
    "morphio_shared. Its target names may have changed in "
    "${OCCUPANCY_MORPHIO_TAG}.")
endif()
