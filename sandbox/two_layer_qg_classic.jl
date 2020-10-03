using QuasiGeostrophy
using LinearAlgebra, Test, FFTW, BenchmarkTools, Plots

# boiler plate definitions
Ωxy = Torus(0,2π) × Torus(0,2π)
Nx = 2^5; Ny = 2^5;
fourier_grid = create_grid((Nx, Ny), Ωxy)
x, y = fourier_grid.grid
kx, ky = fourier_grid.wavenumbers
ψ1 = @. 0 * sin(x) - sin(y) + 0im
ψ2 = @. 0 * sin(x) - sin(y) + 0im
ψ̂1 = copy(ψ1)
ψ̂2 = copy(ψ2)
q̂1 = copy(ψ2)
q̂2 = copy(ψ2)
q1 = copy(ψ2)
q2 = copy(ψ2)
FFTW.set_num_threads(Threads.nthreads())
P = plan_fft(ψ1) # typeof(iP) <: AbstractFFTs.ScaledPlan
iP = plan_ifft(ψ1)
mul!(ψ̂1, P, ψ1)

## Define Parameters
const λ  = 1e-2
const U  = 1e-2 
const ν  = 1e-6
const κ  = 1e-0 
const Δt = 1e-1
const L  = 1.0
const Q  = 1.0

## Define forcing
forcing = @. 0 * sin(x) + Q * sin(y/L) + 0im
f̂ = copy(forcing)
mul!(f̂, P, forcing)

## Define Operators
∂x = FourierOperator(im .* kx)
∂y = FourierOperator(im .* ky)
Δ = ∂x^2 + ∂y^2
Δ⁵ = Δ^5
G = inv( (Δ + -1 /(2*λ^2) )^2 + -(1/(2*λ^2))^2 )

function get_pv!(q̂1, q̂2, ψ̂1, ψ̂2)
    tmp1 = Δ(ψ̂1) + 1 /(2*λ^2) * (ψ̂2 - ψ̂1)
    q̂1 .= tmp1
    tmp2 = Δ(ψ̂2) + 1 /(2*λ^2) * (ψ̂1 - ψ̂2)
    q̂2 .= tmp2
    return nothing
end

get_pv!(q̂1, q̂2, ψ̂1, ψ̂2)

function get_stream_function!(q1, q2, ψ1, ψ2)
    tmp1 = G(Δ(q1) + 1 /(2*λ^2) * (-q2 - q1))
    ψ1 .= tmp1
    tmp2 = G(Δ(q2) + 1 /(2*λ^2) * (-q1 - q2))
    ψ2 .= tmp2
    return nothing
end

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
        q1 = state[1]
        q2 = state[2]
        ψ1 = copy(q1)
        ψ2 = copy(q2)
        get_stream_function!(q1, q2, ψ1, ψ2)
        drag = - 2 .* κ .* (Δ(ψ2))
        q1_rhs = -J(ψ1, q1) - (U/λ^2) .* ∂x(ψ1) - U .* ∂x(q1) + forcing
        q2_rhs = -J(ψ2, q2) + (U/λ^2) .* ∂x(ψ2) + U .* ∂x(q2) - forcing + drag
        return [q1_rhs, q2_rhs]
    end
end


rhs = closure_rhs(∂x, ∂y, Δ, Δ⁵, J, λ, ν, U, κ, f̂)

function closure_lhs(∂x, ∂y, Δ, Δ⁵, J, λ, ν, U, Δt)
    function lhs()
        q1_lhs = 1 + Δt * ν * Δ^4
        return [q1_lhs, q1_lhs]
    end
end

lhs = closure_lhs(∂x, ∂y, Δ, Δ⁵, J, λ, ν, U, Δt)
operators = lhs()
lhs_operator = [inv(operator) for operator in operators]

function imex_step!(state, lhs_operator, rhs, Δt)
    rhs_state = rhs(state)
    current_state = copy(state)
    for i in eachindex(state)
        tmp = lhs_operator[i](current_state[i] .+ Δt .* rhs_state[i])
        @. state[i] =  tmp
    end
end
# @btime imex_step!(state, lhs_operator, rhs, Δt);

##
q1 = @. 1 * sin(x) - cos(y) + 0im
q2 = @. 1 * sin(x) + sin(y) + 0im
fq1 = P * q1
fq2 = P * q2
τ = (q1 - q2) ./ 2
τ̂ = (fq1 - fq2) ./ 2
tmp = P * (400 * sin.(x) .+ 0 * sin.(y))
#fq1 .= good_start1
#fq2 .= good_start2
state = [fq1, fq2]
##
gr(size = (400,400))
for i in 1:2000
    imex_step!(state, lhs_operator, rhs, Δt)
    if (i%100) == 0
        @. τ̂ = (fq1 - fq2) / 2
        mul!(τ, iP, τ̂)
        local p1 = contourf(x[:], y[:], real.(τ)', linewidth = 0, color = :thermometer)
        local p2 = plot(sum(real.(τ), dims = 1)[:] ./ Nx, x[:] )
        display(plot(p1,p2))
        sleep(0.01)
    end 
end
#good_start1 = copy(fq1)
#good_start2 = copy(fq2)
##

@. τ̂ = (fq1 - fq2) / 2
mul!(τ, iP, τ̂)


##

# Check PV and Stream Function Calculation
tmp = @. 1 * sin(2*x) + Q * sin(y/L * 3) + 0im
tmp2 = @. 1 * sin(x) - Q * sin(y/L) + 0im
ψ̂1 = P * tmp
ψ̂2 = P * tmp2
q1 = copy(q̂1)
q2 = copy(q̂2)
ψ1 = copy(ψ̂1)
ψ2 = copy(ψ̂2)

get_pv!(q̂1, q̂2, ψ̂1, ψ̂2)
get_stream_function!(q̂1, q̂2, ψ̂1, ψ̂2)
norm(ψ1 - ψ̂1) 
norm(ψ2 - ψ̂2)
##
i = 1
@. state[i] *= lhs_operator[i].op