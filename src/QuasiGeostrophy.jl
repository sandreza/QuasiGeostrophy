"""
QuasiGeostrophy

# Description

- A series of functions for testing
"""
module QuasiGeostrophy
# Using
using LinearAlgebra
# Spatial Discretization Includes
include("spatial_discretization/nodal_chebyshev.jl")
include("spatial_discretization/fourier_core.jl")
# Spatial Discretization Includes
export cheb
export fourier_nodes, fourier_wavenumbers

# Utility Includes
include("utils/utils.jl")
# Utility Exports
export appropriate_dims, plot

# Symbolic Includes
include("symbolics/abstract_core.jl")

# Grid Includes
include("grid/grid_core.jl")

# Include Hooks to other packages
include("utils/fourier_fields_hooks.jl")
include("utils/fourier_operators_hooks.jl")
include("utils/plot_hooks.jl")
include("utils/convenience_functions.jl")

# Most Important functions
function print_colors()
    for i in 0:255
        printstyled("color(" * string(i) * ")", color = i)
    end
    return nothing
end
export print_colors

end # module
