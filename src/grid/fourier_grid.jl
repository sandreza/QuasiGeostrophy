export FourierGrid
export create_grid

struct FourierGrid{𝒢, 𝒲, 𝒟} <: AbstractGrid
    grid::𝒢
    wavenumbers::𝒲
    domain::𝒟
end

function create_grid(grid_points, Ω::IntervalDomain)
    @assert length(grid_points) == 1
    @assert Ω.periodic
    grid = fourier_nodes(grid_points, a = Ω.a, b = Ω.b)
    wavenumbers = fourier_wavenumbers(grid_points, L = Ω.b - Ω.a)
    return FourierGrid(grid, wavenumbers, Ω)
end

function create_grid(grid_points, Ω::ProductDomain)
    @assert length(grid_points) == length(Ω.domains)
    grid = []
    wavenumbers = []
    @assert check_full_periodicity(Ω)
    for (i,domain) in enumerate(Ω.domains)
        L = domain.b - domain.a
        push!(grid, fourier_nodes(grid_points[i], a = domain.a, b = domain.b))
        push!(wavenumbers, fourier_wavenumbers(grid_points[i], L = L))
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
