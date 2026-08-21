include(FetchContent)

FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG v1.17.0
  GIT_SHALLOW TRUE
  SYSTEM)

set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

add_executable(Tests
               Tests.cpp
               TestAkima.cpp
               TestConcepts.cpp
               TestLagrange.cpp
               TestLinearDeterministic.cpp
               TestPolynomial.cpp)

target_link_libraries(Tests PRIVATE Interpolation::Interpolation gtest_main)
target_compile_options(Tests PRIVATE ${INTERPOLATION_WARNINGS})

include(GoogleTest)
gtest_discover_tests(Tests)
