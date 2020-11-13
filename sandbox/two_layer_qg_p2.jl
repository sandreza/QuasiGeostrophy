using QuasiGeostrophy
using LinearAlgebra, Test, FFTW, BenchmarkTools, Plots
include(pwd() * "/sandbox/qg_utils.jl") # timestepping and wrappers here
# Domain
Ω = S¹(0, 2π) × S¹(0, 2π)
Nx = Ny = 2^5; 
array = Array # to switch between gpu and cpu
fourier_grid = FourierGrid((Nx, Ny), Ω, arraytype = array)
# Fields
fieldnames = ("q¹", "q²", "ψ¹", "ψ²", "f")
q¹, q², ψ¹, ψ², f = create_fields(names = fieldnames, grid = fourier_grid, arraytype = array)
# Operators
∂x, ∂y = create_operators(fourier_grid, arraytype = array)
Δ = ∂x*∂x + ∂y*∂y  # multiplication for GPU purposes
Δ⁴ = -Δ  #* -Δ * -Δ * -Δ  # multiplication for GPU purposes
# Parameters (λ/L)² = 0.01 => λ = 0.1 L 
const λ  = 0.01*2π # Rossby deformation radius, √(gravity * depth) / (coriolis force)
const U  = 1e-4    # background velocity in the top layer, - in bottom
const κ  = 1.0     # bottom drag
const Q  = 0.5     # surface temperature flux

Δtᵐᵃˣ = 1e-1 / Nx     # Largest timestep size
ν = 1e-0/maximum(abs.((Δ⁴.op)))/Δtᵐᵃˣ  # largest stable damping
# initialize stream functions in layers
x, y = fourier_grid.grid
ψ¹( 1.0 * (sin.(x) .* sin.(y) .- 1.0 * sin.(y)) )
ψ²( 1.0 *( (cos.(x) .* sin.(y)) .+ 0.1 * sin.(8*y) .+ 1.0 * sin.(y)  ))
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
wavenums = sqrt.(abs.((dx^2 + dy^2).op))
maxk = maximum(wavenums)
op = cutoff.(wavenums, 1/2 * maxk) #max k is about \sqrt(2) larger than max kx
cutoff_filter = FourierOperator(op, fmd)
ℱ = Operator(nothing, OperatorMetaData(cutoff_filter, "ℱ"))

# more dissipative QG
@pde_system qg_system = (
    ψ¹ = G(Δ(q¹) + 1/(2*λ^2) * (-q²-q¹)),
    ψ² = G(Δ(q²) + 1/(2*λ^2) * (-q¹-q²)),
    ∂t(q¹) = -ν*Δ⁴(q¹) + (ℱ(∂y(ψ¹))*ℱ(∂x(q¹)) - ℱ(∂x(ψ¹))*ℱ(∂y(q¹))) - (U/λ^2) * ∂x(ψ¹) - U * ∂x(q¹) + f,
    ∂t(q²) = -ν*Δ⁴(q²) + (ℱ(∂y(ψ²))*ℱ(∂x(q²)) - ℱ(∂x(ψ²))*ℱ(∂y(q²))) + (U/λ^2) * ∂x(ψ²) + U * ∂x(q²) - f + drag
)

## Timestepping in qg_utils
vort1 = []
vort2 = []
for i in 1:10000
    # adaptive timestepping
    Δt = 0.3 * cfl(ψ¹, ψ², ∂x, ∂y)
    Δt = minimum([Δt, Δtᵐᵃˣ])
    evolve_system(qg_system, Δt)

    if i%100==0
        println("hi at "* string(i))
        
        p1 = plot(q¹.data)
        p2 = plot(q².data)
        display(plot(p1,p2))
        sleep(0.01)
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
        push!( vort1 , real.(P * q¹.data.data))
        push!( vort2 , real.(P * q².data.data))      
    end
end
@save "qg_stuff.jld2" vort1 vort2

plot(q¹.data)
plot(q².data)

spectrum(q¹.data-q².data)
pyplot()
tmp = q¹.data.metadata.transform.backward * (q¹.data.data - q².data.data)
contourf(real.(tmp'), color = :thermometer, linewidth = 0, contourlevlels = 30)
y = sum(real.(tmp'), dims = 2) ./Nx
plot(x, y, label = "horizontal average")
plot!(x,sin.(x) * maximum(y), label = "sin")