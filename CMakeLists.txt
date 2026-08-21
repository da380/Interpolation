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
# The library
# ---------------------------------------------------------------------------
add_library(Interpolation INTERFACE)
add_library(Interpolation::Interpolation ALIAS Interpolation)

target_compile_features(Interpolation INTERFACE cxx_std_23)

target_include_directories(Interpolation INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

# ---------------------------------------------------------------------------
# Install and export, so find_package(Interpolation) works from an install tree
# ---------------------------------------------------------------------------
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
