using QuasiGeostrophy
using LinearAlgebra, Test, FFTW, BenchmarkTools, Plots, JLD2, CUDA
include(pwd() * "/sandbox/qg_utils.jl") # timestepping and wrappers here
# Domain
Ω = S¹(0, 2π) × S¹(0, 2π)
Nx = Ny = 2^8; 
array = CuArray # to switch between gpu and cpu
fourier_grid = FourierGrid((Nx, Ny), Ω, arraytype = array)
# Fields
fieldnames = ("q¹", "q²", "ψ¹", "ψ²", "f")
q¹, q², ψ¹, ψ², f = create_fields(names = fieldnames, grid = fourier_grid, arraytype = array)
# Operators
∂x, ∂y = create_operators(fourier_grid, arraytype = array)
Δ = ∂x*∂x + ∂y*∂y  # multiplication for GPU purposes
Δ⁴ = -Δ * -Δ #* -Δ * -Δ  # multiplication for GPU purposes
# Parameters
const λ  = 0.01 * 2π   
const U  = 0e-0 
const κ  = 1.0
Δt = 1e-1 / Nx
const L  = 1.0
const Q  = 0.5  
ν = 1e-0/maximum(abs.((Δ⁴.op)))/Δt # largest stable damping
# initialize stream functions in layers
x, y = fourier_grid.grid
ψ¹( 0.1 * (sin.(x) .* sin.(y) .- 1.0 * sin.(y)) )
ψ²( 0.1 *( (cos.(x) .* sin.(y)) .+ 0.1 * sin.(8*y) .+ 1.0 * sin.(y)  ))
# initialize forcing
f(Q * sin.(y))
# stream function to potential vorticity
q¹.data .= (Δ(ψ¹) + 1 /(2*λ^2) * (ψ² - ψ¹)).data
q².data .= (Δ(ψ²) + 1 /(2*λ^2) * (ψ¹ - ψ²)).data
# potential vorticity to stream function
G = inv( (Δ + -1 /(2*λ^2) ) * (Δ + -1 /(2*λ^2) ) + -(1/(2*λ^2))^2 )
ψ¹.data .= G(Δ(q¹) + 1 /(2*λ^2) * (-q² - q¹)).data
ψ².data .= G(Δ(q²) + 1 /(2*λ^2) * (-q¹ - q²)).data

# Impero
@wrapper q¹=q¹ ψ¹=ψ¹ q²=q² ψ²=ψ² f=f
Δ = Operator(Δ); ∂x = Operator(∂x); ∂y = Operator(∂y)
Δ⁴ = Operator(Δ⁴);
G = Operator(nothing, OperatorMetaData(G, "G"))
∂t = Operator(nothing, OperatorMetaData(nothing, "∂t"))
drag = (-2 * κ) * Δ(ψ²)

# create cutoff filter
fmd = FourierOperatorMetaData("Filter")
dx, dy = create_operators(fourier_grid, arraytype = array)
wavenums = sqrt.(abs.((dx*dx + dy*dy).op))
maxk = maximum(wavenums)
op = cutoff.(wavenums, 1/3 * maxk) #max k is about \sqrt(2) larger than max kx
cutoff_filter = FourierOperator(op, fmd)
ℱ = Operator(nothing, OperatorMetaData(cutoff_filter, "ℱ"))

# more dissipative QG
@pde_system qg_system = (
    ψ¹ = G(Δ(q¹) + 1/(2*λ^2) * (-q²-q¹)),
    ψ² = G(Δ(q²) + 1/(2*λ^2) * (-q¹-q²)),
    ∂t(q¹) = -ν*Δ⁴(q¹) + (ℱ(∂y(ψ¹))*ℱ(∂x(q¹)) - ℱ(∂x(ψ¹))*ℱ(∂y(q¹))) - (U/λ^2) * ∂x(ψ¹) - U * ∂x(q¹) + f,
    ∂t(q²) = -ν*Δ⁴(q²) + (ℱ(∂y(ψ¹))*ℱ(∂x(q¹)) - ℱ(∂x(ψ¹))*ℱ(∂y(q¹))) + (U/λ^2) * ∂x(ψ²) + U * ∂x(q²) - f + drag
)

## Timestepping in qg_utils
vort1 = []
vort2 = []
for i in 1:100000
    # adaptive timestepping
    Δt = 0.1 * cfl(ψ¹, ψ², ∂x, ∂y)
    Δt = minimum([Δt, 1e-1 / Nx])
    evolve_system(qg_system, Δt, filter = cutoff_filter)
    if (i%1000==0) && (i > 20000)
        #=
        p1 = plot(q¹.data)
        p2 = plot(q².data)
        display(plot(p1,p2))
        sleep(0.01)
        =#
        println("hi at "* string(i))
        println("The norm is ")
        println(norm(compute(∂y(ψ¹)).data)/length(ψ¹.data.data))
        println("The norm is ")
        println(norm(compute(∂y(ψ²)).data)/length(ψ².data.data))
        println("The Δt is ")
        println(Δt)
        println("The CFL is ")
        println(Δt/cfl(ψ¹, ψ², ∂x, ∂y) )
        P = q¹.data.metadata.transform.backward
        push!( vort1, Array(real.(P * q¹.data.data)))
        push!( vort2, Array(real.(P * q².data.data)))      
    end
end
@save "qg_stuff.jld2" vort1 vort2