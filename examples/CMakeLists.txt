# The examples are numbered in the order they are meant to be read: each one
# assumes what came before it. They print to stdout rather than writing data
# files, so running one tells you immediately whether it did the right thing.

set(INTERPOLATION_EXAMPLES
    01-linear-interpolation
    02-cubic-splines
    03-boundary-conditions
    04-akima-and-lagrange
    05-polynomials
    06-borrowing-and-owning
    07-function-algebra
    08-derivatives-and-integrals
    09-piecewise-layers
    10-two-dimensions)

foreach(example IN LISTS INTERPOLATION_EXAMPLES)
  add_executable(${example} ${example}.cpp)
  target_link_libraries(${example} PRIVATE Interpolation::Interpolation)
  target_compile_options(${example} PRIVATE ${INTERPOLATION_WARNINGS})
endforeach()

# One target to build the lot.
add_custom_target(examples DEPENDS ${INTERPOLATION_EXAMPLES})
