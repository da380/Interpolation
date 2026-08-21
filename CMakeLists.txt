cmake_minimum_required(VERSION 3.24)

project(Interpolation
        VERSION 0.1.0
        DESCRIPTION "Header-only C++23 library for one-dimensional interpolation and polynomial arithmetic"
        HOMEPAGE_URL "https://github.com/da380/Interpolation"
        LANGUAGES CXX)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
  message(FATAL_ERROR
    "In-source builds are not allowed. Configure into a separate build "
    "directory, for example: cmake -S . -B build. You may need to remove "
    "CMakeCache.txt and CMakeFiles/ first.")
endif()

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# ---------------------------------------------------------------------------
# Dependencies
#
# Eigen is transitional: cubic splines still use SimplicialLDLT on a system
# that is only tridiagonal, and AkimaSpline includes three Eigen headers it
# never uses. Phase 2 replaces the solve with a Thomas algorithm and removes
# the dependency entirely.
# ---------------------------------------------------------------------------
# No version is requested, deliberately. Eigen 3.4 installs a SameMajorVersion
# version file that rejects any range spanning major versions, while Eigen 5
# rejects a bare "3.4" as meaning 3.4.x only, so no single versioned request
# satisfies both. Everything this library uses has been present since 3.3, and
# the dependency goes away in phase 2 regardless.
#
# The package registry is excluded deliberately: it otherwise resolves Eigen3
# from an unrelated sibling project's build tree, which makes local builds
# depend on directories outside this repository and diverge from CI.
find_package(Eigen3 QUIET NO_MODULE
             NO_CMAKE_PACKAGE_REGISTRY
             NO_CMAKE_SYSTEM_PACKAGE_REGISTRY)

if(NOT Eigen3_FOUND)
  include(FetchContent)
  FetchContent_Declare(
    Eigen3
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 5.0.1
    GIT_SHALLOW TRUE
    SYSTEM)
  set(EIGEN_BUILD_TESTING OFF CACHE BOOL "" FORCE)
  set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
  set(EIGEN_BUILD_DEMOS OFF CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(Eigen3)
endif()

# ---------------------------------------------------------------------------
# The library
# ---------------------------------------------------------------------------
add_library(Interpolation INTERFACE)
add_library(Interpolation::Interpolation ALIAS Interpolation)

target_compile_features(Interpolation INTERFACE cxx_std_23)

target_link_libraries(Interpolation INTERFACE Eigen3::Eigen)

target_include_directories(Interpolation INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

# ---------------------------------------------------------------------------
# Install and export, so find_package(Interpolation) works from an install tree
#
# A FetchContent-supplied Eigen is not part of our export set, so
# install(EXPORT) cannot reference it and the package is not installable in
# that configuration. Installing therefore requires a real, findable Eigen3.
# Phase 2 removes Eigen and this condition with it.
# ---------------------------------------------------------------------------
if(Eigen3_FOUND)
  install(TARGETS Interpolation EXPORT InterpolationTargets)

  install(DIRECTORY include/Interpolation
          DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
          FILES_MATCHING PATTERN "*.hpp")

  install(EXPORT InterpolationTargets
          FILE InterpolationTargets.cmake
          NAMESPACE Interpolation::
          DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Interpolation)

  write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/InterpolationConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
    ARCH_INDEPENDENT)

  configure_package_config_file(
    cmake/InterpolationConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/InterpolationConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Interpolation)

  install(FILES
          ${CMAKE_CURRENT_BINARY_DIR}/InterpolationConfig.cmake
          ${CMAKE_CURRENT_BINARY_DIR}/InterpolationConfigVersion.cmake
          DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Interpolation)
else()
  message(STATUS
    "Eigen3 was fetched rather than found, so it is not in this project's "
    "export set and install rules are disabled. Point CMAKE_PREFIX_PATH at an "
    "installed Eigen3 to build an installable package.")
endif()

# ---------------------------------------------------------------------------
# Developer targets, only when this is the top-level project
# ---------------------------------------------------------------------------
if(PROJECT_IS_TOP_LEVEL)
  option(INTERPOLATION_BUILD_EXAMPLES "Build the example programs" ON)
  option(INTERPOLATION_BUILD_TESTS    "Build the test suite" ON)
  option(INTERPOLATION_BUILD_DOCS     "Build the Doxygen API documentation" OFF)

  option(INTERPOLATION_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)

  # Warnings are enabled on our own targets rather than globally, so that
  # FetchContent dependencies are not held to them.
  if(MSVC)
    set(INTERPOLATION_WARNINGS /W4)
    if(INTERPOLATION_WARNINGS_AS_ERRORS)
      list(APPEND INTERPOLATION_WARNINGS /WX)
    endif()
  else()
    set(INTERPOLATION_WARNINGS -Wall -Wextra -Wpedantic)
    if(INTERPOLATION_WARNINGS_AS_ERRORS)
      list(APPEND INTERPOLATION_WARNINGS -Werror)
    endif()
  endif()

  if(INTERPOLATION_BUILD_EXAMPLES)
    add_subdirectory(examples)
  endif()

  if(INTERPOLATION_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
  endif()

  if(INTERPOLATION_BUILD_DOCS)
    find_package(Doxygen REQUIRED)
    set(INTERPOLATION_DOXYGEN_CONFIG ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
    configure_file(docs/Doxyfile.in ${INTERPOLATION_DOXYGEN_CONFIG} @ONLY)
    add_custom_target(InterpolationDocs
      COMMAND ${DOXYGEN_EXECUTABLE} ${INTERPOLATION_DOXYGEN_CONFIG}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      COMMENT "Generating Interpolation API documentation"
      VERBATIM)
  endif()
endif()
