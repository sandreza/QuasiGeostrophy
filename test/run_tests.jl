include(pwd() * "/test/ode_solution.jl")

# Impero Tests
include(pwd() * "/test/symbolic_tests.jl")
include(pwd() * "/test/symbolic_utils_tests.jl")
include(pwd() * "/test/test_equation.jl")

# Fourier Objects Tests
include(pwd() * "/test/fourier_tests.jl")
include(pwd() * "/test/fourier_field_tests.jl")

# Fourier Impero Interaction Test
include(pwd() * "/test/impero_fourier_test.jl")

# Timestepping Test
include(pwd() * "/test/time_step_test.jl")