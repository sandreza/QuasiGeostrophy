using QuasiGeostrophy
using LinearAlgebra, Test, FFTW, BenchmarkTools, Plots
include(pwd() * "/sandbox/qg_utils.jl") # timestepping and wrappers here
# Domain
Ω = S¹(0, 2π) × S¹(0, 2π)
Nx = Ny = 2^6; 
array = Array # to switch between gpu and cpu
fourier_grid = FourierGrid((Nx, Ny), Ω, arraytype = array)
# Fields
fieldnames = ("q¹", "q²", "ψ¹", "ψ²", "f")
q¹, q², ψ¹, ψ², f = create_fields(names = fieldnames, grid = fourier_grid, arraytype = array)
# Operators
∂x, ∂y = create_operators(fourier_grid, arraytype = array)
Δ = ∂x*∂x + ∂y*∂y  # multiplication for GPU purposes
Δ⁴ = Δ * Δ * Δ * Δ # multiplication for GPU purposes
# Parameters
const λ  = 1e-1
const U  = 1e-0 
const ν  = 3e-6 / λ^2 / Nx
const κ  = 1e-0 
const Δt = 1e-1 / Nx
const L  = 2π
const Q  = 1.0
# filter
filter = inv(1 + ν*Δ⁴)
# initialize stream functions in layers
x, y = fourier_grid.grid
ψ¹(0 .* sin.(x) .+ sin.(y))
ψ²(cos.(x) .* sin.(y))
# initialize forcing
f(Q * sin.(y .* 2π/L) .+ 0.0 .* x)
# stream function to potential vorticity
q¹.data .= (Δ(ψ¹) + 1 /(2*λ^2) * (ψ² - ψ¹)).data
q².data .= (Δ(ψ²) + 1 /(2*λ^2) * (ψ¹ - ψ²)).data
# potential vorticity to stream function
G = inv( (Δ + -1 /(2*λ^2) ) * (Δ + -1 /(2*λ^2) ) + -(1/(2*λ^2))^2 )
ψ¹.data .= G(Δ(q¹) + 1 /(2*λ^2) * (-q² - q¹)).data
ψ².data .= G(Δ(q²) + 1 /(2*λ^2) * (-q¹-q²)).data

# Impero
@wrapper q¹=q¹ ψ¹=ψ¹ q²=q² ψ²=ψ² f=f
Δ = Operator(Δ); ∂x = Operator(∂x); ∂y = Operator(∂y)
G = Operator(nothing, OperatorMetaData(G, "G"))
∂t = Operator(nothing, OperatorMetaData(nothing, "∂t"))
drag = (-2 * κ) * Δ(ψ²)

@pde_system qg_system = (
    ψ¹ = G(Δ(q¹) + 1/(2*λ^2) * (-q²-q¹)),
    ψ² = G(Δ(q²) + 1/(2*λ^2) * (-q¹-q²)),
    ∂t(q¹) = (∂y(ψ¹)*∂x(q¹) - ∂x(ψ¹)*∂y(q¹)) - (U/λ^2) * ∂x(ψ¹) - U * ∂x(q¹) + f,
    ∂t(q²) = (∂y(ψ²)*∂x(q²) - ∂x(ψ²)*∂y(q²)) + (U/λ^2) * ∂x(ψ²) + U * ∂x(q²) - f + drag
)

## Timestepping in qg_utils
for i in 1:10000
    evolve_system(qg_system, Δt, filter = filter)
    if i%1000==0
        p1 = plot(q¹.data)
        display(p1)
        sleep(0.01)
        println("hi at "* string(i))
    end
end

plot(q¹.data)

spectrum(q¹.data+q².data)

tmp = q¹.data.metadata.transform.backward * (q¹.data.data - q².data.data)
contourf(real.(tmp'), color = :thermometer)
y = sum(real.(tmp'), dims = 2) ./Nx
plot(x, y, label = "horizontal average")
plot!(x,- sin.(x) * maximum(y), label = "sin")