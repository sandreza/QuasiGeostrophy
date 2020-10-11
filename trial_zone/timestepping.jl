using QuasiGeostrophy, FFTW, Test, BenchmarkTools
using LinearAlgebra
import QuasiGeostrophy: to_expr
include(pwd() * "/test/test_utils.jl")
const tol = 1e1

## Define 1D Test
Ω = S¹(0, 2π) 
Nx = 2^7; 
fourier_grid = FourierGrid(Nx, Ω)
fieldnames = ("u", "σ", "v")
u, σ, v = create_fields(names = fieldnames, grid = fourier_grid)
∂x = create_operators(fourier_grid)
# initialize FourierField data
x = fourier_grid.grid[1]
u(sin.(x) .+ 1); σ(cos.(x));
v( sin.(x));
# Wrap Around Impero Object
@wrapper u=u v=v σ=σ
Δ  = Operator(∂x^2)
∂x = Operator(∂x)
κ = 0.0
∂t = Operator(nothing, OperatorMetaData(nothing, "∂t"))

@pde_system pde_system = [
    ∂t(u) = -∂x(u) + κ*Δ(u),
]
to_expr(eq::Equation) = Expr(:(=), to_expr(eq.lhs), to_expr(eq.rhs))
plot(to_expr(pde_system[1]))
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

for i in 1:1000
    evolve_pde(pde_system, 0.01, RK2())
    if i%10==0
        p1 = plot(u.data)
        display(p1)
        sleep(0.01)
    end
end