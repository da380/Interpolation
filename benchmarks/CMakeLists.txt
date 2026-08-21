add_executable(InterpolationBenchmarks Benchmarks.cpp)
target_link_libraries(InterpolationBenchmarks
                      PRIVATE Interpolation::Interpolation)
target_compile_options(InterpolationBenchmarks PRIVATE ${INTERPOLATION_WARNINGS})
