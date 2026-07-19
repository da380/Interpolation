include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG main
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

add_executable(Tests
               Tests.cpp
               TestAkima.cpp
               TestConcepts.cpp
               TestLagrange.cpp
               TestLinearDeterministic.cpp
               TestPolynomial.cpp)
target_link_libraries(Tests PRIVATE Interpolation gtest_main)
include(GoogleTest)
gtest_discover_tests(Tests)
