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
    include("spatial_discretization/fourier.jl")
    # Spatial Discretization Includes
    export cheb
    export fourier_nodes, fourier_wavenumbers
    
    # Utility Includes
    include("utils/utils.jl")
    # Utility Exports
    export appropriate_dims
end # module
