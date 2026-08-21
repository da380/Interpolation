cmake_minimum_required (VERSION 3.5)



project(Interpolation)


if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt. ")
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED true)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# add in Eigen3
include(FetchContent)
FetchContent_Declare(
  Eigen3
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG master
)
FetchContent_MakeAvailable(Eigen3)

# Set up the library
add_library (Interpolation INTERFACE)
target_link_libraries(Interpolation INTERFACE Eigen3::Eigen)
target_include_directories (Interpolation INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${INCLUDE_INSTALL_DIR}>
)


# optionally add in the examples and tests
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    # We're in the root, define additional targets for developers.
    option(MY_PROJECT_BUILD_EXAMPLES   "whether or not examples should be built" ON)
    option(MY_PROJECT_BUILD_TESTS      "whether or not tests should be built" ON)
    option(INTERPOLATION_BUILD_DOCS    "whether or not API documentation should be built" OFF)

    if(MY_PROJECT_BUILD_EXAMPLES)
        add_subdirectory(examples)
    endif()
    if(MY_PROJECT_BUILD_TESTS)
        enable_testing()
        add_subdirectory(tests)
    endif()
    if(INTERPOLATION_BUILD_DOCS)
        find_package(Doxygen REQUIRED)
        set(INTERPOLATION_DOXYGEN_CONFIG
            ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
        configure_file(docs/Doxyfile.in
                       ${INTERPOLATION_DOXYGEN_CONFIG}
                       @ONLY)
        add_custom_target(InterpolationDocs
            COMMAND ${DOXYGEN_EXECUTABLE} ${INTERPOLATION_DOXYGEN_CONFIG}
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            COMMENT "Generating Interp API documentation"
            VERBATIM)
    endif()
endif()

