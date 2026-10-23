add_executable(interpolation_benchmarks Benchmarks.cpp)
set_target_properties(interpolation_benchmarks PROPERTIES
                      OUTPUT_NAME InterpolationBenchmarks)
target_link_libraries(interpolation_benchmarks
                      PRIVATE Interpolation::Interpolation)
target_compile_options(interpolation_benchmarks PRIVATE ${INTERPOLATION_WARNINGS})
