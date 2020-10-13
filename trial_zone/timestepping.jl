using QuasiGeostrophy, FFTW, Test, BenchmarkTools
using LinearAlgebra
using Plots, GraphRecipes
using SymbolicUtils
import SymbolicUtils: Chain, Postwalk
include(pwd() * "/test/test_utils.jl")
const tol = 1e1

## Define 1D Test
Ω = S¹(0, 2π) 
Nx = 2^9; 
fourier_grid = FourierGrid(Nx, Ω)
fieldnames = ("u", "σ", "v")
u, σ, v = create_fields(names = fieldnames, grid = fourier_grid)
∂x = create_operators(fourier_grid)
# initialize FourierField data
x = fourier_grid.grid[1]
u(sin.(x) .+ 1); σ(cos.(x));
v( sin.(x));
starting_energy = (u*u).data[1]
# Wrap Around Impero Objects
@wrapper u=u v=v σ=σ
Δ  = Operator(∂x^2)
∂x = Operator(∂x)
κ = 10.0/Nx

# This is just syntax for now
∂t = Operator(nothing, OperatorMetaData(nothing, "∂t"))

@pde_system pde_system = [
    ∂t(u) = -∂x(u^2) + κ*Δ(u),
]

pde_plot = plot(pde_system[1]);
abstract type TimeSteppingMethod end
struct  RK1 <: TimeSteppingMethod end #include Δt as a part of RK1
struct  RK2 <: TimeSteppingMethod end

function evolve_pde(pde_system, Δt, ::RK1)
    data = pde_system[1].lhs.operand.data.data
    rhs = evaluate(pde_system[1].rhs)
    data .+= Δt .* rhs
    return nothing
end

function evolve_pde(pde_system, Δt, ::RK2)
    data = pde_system[1].lhs.operand.data.data
    olddata = copy(data)
    stage1 = evaluate(pde_system[1].rhs)
    data .+= Δt .* stage1
    stage2 = evaluate(pde_system[1].rhs)
    data .= olddata .+ Δt/2 .* (stage1 + stage2)
    return nothing
end

const Δt = 0.1 / Nx  
T = 1
steps = round(Int, T / Δt) 
plotsteps = round(Int, steps / 10)
for i in 1:steps
    evolve_pde(pde_system, Δt, RK2())
    if i%plotsteps==0
        p1 = plot(u.data, ylims = (0,2))
        display(p1)
        sleep(0.01)
    end
end

ending_energy = evaluate(u * u)[1]
energy_loss = (starting_energy - ending_energy)/starting_energy * 100
