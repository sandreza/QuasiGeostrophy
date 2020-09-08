using QuasiGeostrophy, LinearAlgebra, Test, FFTW, BenchmarkTools, Plots

# boiler plate definitions
Ωxy = Torus(0,2π) × Torus(0,2π)
Nx = 2^7; Ny = 2^7;
fourier_grid = create_grid((Nx, Ny), Ωxy)
x, y = fourier_grid.grid
kx, ky = fourier_grid.wavenumbers
ψ = @. 0 * sin(x) - sin(y) + 0im
τ = @. 0 * sin(x) - sin(y) + 0im
ψ̂ = copy(ψ)
τ̂ = copy(τ)
FFTW.set_num_threads(Threads.nthreads())
P = plan_fft(ψ) # typeof(iP) <: AbstractFFTs.ScaledPlan
iP = plan_ifft(ψ)
mul!(ψ̂, P, ψ)

## Define Parameters
const λ  = 0.1 
const U  = 1.0 
const ν  = 1.0 
const κ  = 1.0 
const Δt = 1e-2 
const L  = 1.0
const Q  = 100.0

## Define forcing
forcing = @. 0 * sin(x) + Q * sin(y/L) + 0im
f̂ = copy(forcing)
mul!(f̂, P, forcing)

## Define Operators
∂x = FourierDerivative(im .* kx)
∂y = FourierDerivative(im .* ky)
Δ = ∂x^2 + ∂y^2
Δ⁵ = Δ^5

function closure_J(∂x, ∂y, P, iP)
    function J(a,b)
        term_1 = filter_convolve(∂x(a) , ∂y(b), P, iP)
        term_2 = filter_convolve(∂x(b) , ∂y(a), P, iP)
        return  term_1 - term_2
    end
end

J = closure_J(∂x, ∂y, P, iP)

function closure_rhs(∂x, ∂y, Δ, Δ⁵, J, λ, ν, U, κ, forcing)
    function rhs(state)
        ψ = state[1]
        τ = state[2]
        Δψ = Δ(ψ)
        Δτ = Δ(τ)
        λ⁻²τ = τ ./ (λ^2)
        λ⁻²ψ = ψ ./ (λ^2)
        drag = - 2 .* κ .* (Δψ - Δτ) # check this
        ψ_rhs = -J(ψ, Δψ) - J(τ, Δτ) - U .* ∂x(Δτ) + drag .* 0.5 
        τ_rhs = -J(ψ, Δτ - λ⁻²τ) - J(τ, Δψ) - U .* ∂x(Δψ + λ⁻²ψ) - drag .* 0.5 + forcing 
        return [ψ_rhs, τ_rhs]
    end
end

rhs = closure_rhs(∂x, ∂y, Δ, Δ⁵, J, λ, ν, U, κ, f̂)

function closure_lhs(∂x, ∂y, Δ, Δ⁵, J, λ, ν, U, Δt)
    function lhs(ψ, τ)
        ψ_lhs = Δ + Δt * ν * Δ⁵
        τ_lhs = ( Δ + -(1/λ^2) ) + Δt * ( ν * Δ⁵ + -(ν/λ^2) * (Δ^4) )
        return [ψ_lhs, τ_lhs]
    end
end

lhs = closure_lhs(∂x, ∂y, Δ, Δ⁵, J, λ, ν, U, Δt)
operators = lhs(ψ,τ)
lhs_operator = [inv(operator) for operator in operators]

function imex_step!(state, lhs_operator, rhs, Δt)
    rhs_state = rhs(state)
    for i in eachindex(state)
        @. state[i] =  lhs_operator[i].op * (state[i] + Δt * rhs_state[i])
    end
end

##
ψ = @. 0 * sin(x) + sin(y) + 0im
τ = @. 0 * sin(x) + sin(y) + 0im
ψ̂ = copy(ψ)
τ̂ = copy(τ)
FFTW.set_num_threads(Threads.nthreads())
P = plan_fft(ψ) # typeof(iP) <: AbstractFFTs.ScaledPlan
iP = plan_ifft(ψ)
mul!(ψ̂, P, ψ)
##
state = [ψ̂, τ̂];
gr(size = (200,200))
##
for i in 1:20000
    imex_step!(state, lhs_operator, rhs, Δt)
    if (i%100) == 0
        mul!(τ, iP, τ̂)
        local p1 = contourf(x[:],y[:], real.(τ)')
        display(p1)
        sleep(0.01)
    end 
end
##
p1 = contourf(x[:], y[:], real.(τ)')
display(p1)
##
i = 1
lhs_operator[i].op .* (state[i] .+ Δt .* rhs(state)[i])