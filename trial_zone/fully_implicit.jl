using QuasiGeostrophy, FFTW, Test, BenchmarkTools
using LinearAlgebra
using Plots, GraphRecipes
using SymbolicUtils
import SymbolicUtils: Chain, Postwalk
include(pwd() * "/test/test_utils.jl")
const tol = 1e1
gr(size = (300,300))
# Define 1D Test
Ω = S¹(0, 2π) 
Nx = Ny = Nz = 2^3; 
fourier_grid = FourierGrid(Nx, Ω)
fieldnames = ("u", "σ", "v", "δu")
u, σ, v, δu = create_fields(names = fieldnames, grid = fourier_grid)
∂x = create_operators(fourier_grid)
Δ = ∂x^2
# initialize FourierField data
x = fourier_grid.grid[1]
Δx = x[2] - x[1]
u(sin.(x) .+ 0.5); σ(cos.(x));
v(sin.(x));
starting_energy = (u*u).data[1]

# initialize parameters
const Δt = 0.1 / Nx  
const ν = 0.0/(2π)^2 * Δx^2 / Δt

# Forward-Euler
# uⁿ⁺¹ = uⁿ + Δt u̇ⁿ
T = 1.5
steps = round(Int, T / Δt) 
plotsteps = round(Int, steps / 10)
for i in 1:steps
    u̇ = - ∂x( 0.5 * u * u) + ν * Δ(u)
    u.data .+= Δt * u̇.data
end
##

# Backward-Euler
# uⁿ⁺¹ = uⁿ + Δt u̇ⁿ⁺¹
# solve uⁿ⁺¹ - Δt *( ∂x( 0.5 * uⁿ⁺¹ * uⁿ⁺¹) - νΔuⁿ⁺¹ ) = uⁿ
# Jacobian
function Jacobian(u, δu, ν, Δt)
    return δu - Δt *( ∂x(u * δu) - ν*Δ(δu) )
end
# Residual
function Residual(uⁿ⁺¹, uⁿ, ν, Δt)
    u̇ⁿ⁺¹ = - ∂x( 0.5 * uⁿ⁺¹ * uⁿ⁺¹) + ν * Δ(uⁿ⁺¹)
    return uⁿ-uⁿ⁺¹ + Δt * u̇ⁿ⁺¹ 
end
# Create Jacobian Matrix for "Exact" Jacobian
function create_matrix(J, δu)
    n = length(δu.data)
    A = zeros(eltype(δu.data),(n,n))
    δu.data .= 0
    for i in 1:n
        for j in [1.0, 1.0im]
            δu.data[i] = j
            jA = J(δu).data
            A[i,:] .+= jA
            δu.data[i] = 0.0
        end
    end
    return A
end

##
u.data .= 0.0
v.data .= 0.0
δu.data .= 0.0
u(sin.(x) .+ 0.5); σ(cos.(x));
v(sin.(x));

const Δt = 1.0 / Nx  
const ν = 0.0/(2π)^2 * Δx^2 / Δt

# Forward-Euler
# uⁿ⁺¹ = uⁿ + Δt u̇ⁿ
T = 1.5
steps = round(Int, T / Δt) 
plotsteps = round(Int, steps / 10)

for i in 1:steps
    v.data .= u.data
    # initial guess for δu
    #δu.data .= Δt * (- ∂x( 0.5 * u * u) + ν * Δ(u)).data
    δu.data .= 0.0
    u.data .+= δu.data
    # newton iteration
    newtonN = 1
    for j in 1:newtonN
        residual = Residual(u, v, ν, Δt)
        #println("residual at ", i, " is ", norm(residual))
        J(δu) = Jacobian(u.data, δu, ν, Δt)
        A = create_matrix(J, δu)
        δu.data .= A \ residual.data
        u.data .+= δu.data
    end
    residual = Residual(u, v, ν, Δt)
    println("residual at ", i, " now is ", norm(residual))
end

