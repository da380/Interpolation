foreach(example
        LinearExample
        CubicSplineExample
        LagrangeExample
        testmult)
  add_executable(${example} ${example}.cpp)
  target_link_libraries(${example} PRIVATE Interpolation::Interpolation)
  target_compile_options(${example} PRIVATE ${INTERPOLATION_WARNINGS})
endforeach()
