export FourierGrid
export create_grid, getindex

struct FourierGrid{𝒢, 𝒲, 𝒟} <: AbstractGrid
    grid::𝒢
    wavenumbers::𝒲
    domain::𝒟
end

import Base: getindex
getindex(f::FourierGrid, i) = f.grid[i]

function create_grid(grid_points, Ω::IntervalDomain) # perhaps change to match the product domain
    @assert length(grid_points) == 1
    @assert Ω.periodic
    grid = fourier_nodes(grid_points, a = Ω.a, b = Ω.b)
    wavenumbers = fourier_wavenumbers(grid_points, L = Ω.b - Ω.a)
    return FourierGrid(grid, wavenumbers, Ω)
end

"""
create_grid(grid_points, Ω::ProductDomain)
# Description
create a numerical grid with grid_points resolution in the domain Ω
# Arguments
- `grid_points`: tuple | a tuple of ints in each direction for product domain
- `Ω`: ProductDomain   | a product domain object
# Return
A Fourier Grid object
"""
function create_grid(grid_points, Ω::ProductDomain) # change to be fully general later
    @assert length(grid_points) == length(Ω.domains)
    grid = []
    wavenumbers = []
    @assert check_full_periodicity(Ω)
    for (i,domain) in enumerate(Ω.domains)
        L = domain.b - domain.a
        reshape_dims = appropriate_dims(length(Ω.domains), i, grid_points[i])
        push!(grid, reshape(fourier_nodes(grid_points[i], a = domain.a, b = domain.b), reshape_dims))
        push!(wavenumbers, reshape(fourier_wavenumbers(grid_points[i], L = L), reshape_dims))
    end
    return FourierGrid(grid, wavenumbers, Ω)
end

function Base.show(io::IO, F::FourierGrid) 
    print("domain:", F.domain, " gridpoints:")
    if typeof(F.grid[1]) <: AbstractArray
        for (i, grid) in enumerate(F.grid)
            print(length(grid))
            if i != length(F.grid)
                printstyled(io, "×")
            end
        end
    else
        print(length(F.grid))
    end
 end
