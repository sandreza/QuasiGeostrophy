using QuasiGeostrophy, CUDA, BenchmarkTools, FFTW

Ω = S¹(0, 2π) × S¹(0, 2π)
Nx = Ny = 2^10; 
arraytype = CuArray
fourier_grid = FourierGrid((Nx, Ny), Ω, arraytype = arraytype )
fieldnames = ("u", "v", "ψ")
u, v, ψ = create_fields(names = fieldnames, grid = fourier_grid, arraytype = arraytype )
∂x, ∂y = create_operators(fourier_grid, arraytype = arraytype )
# initialize FourierField data
x, y = fourier_grid.grid
ψ(sin.(x) .* sin.(y))
u.data .= ∂y(ψ).data
v.data .= -∂x(ψ).data
Δ = FourierOperator((∂x * ∂x + ∂y * ∂y).op, FourierOperatorMetaData("Δ"))
Δ⁻¹ = inv(Δ) 
Re = 0.1
# 12 fourier transforms (3 for each field multiplications)
function nsrhs(u,v, ∂x, ∂y, Δ, Δ⁻¹, Re)
    rhsu = -(∂x(u*u) + ∂y(u*v)) + (1/Re)*Δ(u)
    rhsv = -(∂x(u*v) + ∂y(v*v)) + (1/Re)*Δ(v)
    p = Δ⁻¹(∂x(rhsu) + ∂y(rhsv))
    rhsu -= ∂x(p)
    rhsv -= ∂y(p)
    return rhsu, rhsv
end
@btime nsrhs(u,v, ∂x, ∂y, Δ, Δ⁻¹, Re);