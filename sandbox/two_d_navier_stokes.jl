using QuasiGeostrophy, LinearAlgebra, Test, FFTW, BenchmarkTools, Plots

# boiler plate definitions
Ωxy = Torus(0,2π) × Torus(0,2π)
Nx = 2^6; Ny = 2^6;
fourier_grid = create_grid((Nx, Ny), Ωxy)
x, y = fourier_grid.grid
kx, ky = fourier_grid.wavenumbers
u = @. 1 * sin(x) * sin(y) + 0im
v = @. 1 * cos(x) * cos(y) + 0im
û = copy(u)
v̂ = copy(v)


FFTW.set_num_threads(Threads.nthreads())
P = plan_fft(u) # typeof(iP) <: AbstractFFTs.ScaledPlan
iP = plan_ifft(u)

uv = P * (u .* v)
norm(uv - convolve(P * u, P * v, P, iP)) / norm(uv)

## Define Parameter
const Re  = 100.0 
const Δt = 1.0 / Re
const G  = 4.0

## Define forcing
forcing = @. G * sin(5 * x) * sin(5 * y) + 0im
f̂ = copy(forcing)
mul!(f̂, P, forcing)

## Define Operators
∂x = FourierDerivative(im .* kx)
∂y = FourierDerivative(im .* ky)
Δ = ∂x^2 + ∂y^2
Δ⁻¹ = inv(Δ)


function get_pressure(u_rhs, v_rhs, ∂x, ∂y, G)
    return G(∂x(u_rhs) + ∂y(v_rhs))
end

function closure_rhs(∂x, ∂y, Δ⁻¹, P, iP, forcing)
    function rhs(state)
        u = state[1]
        v = state[2]
        uu = convolve(u, u, P, iP)
        uv = convolve(u, v, P, iP)
        vv = convolve(v, v, P, iP)
        u_rhs = - ∂x(uu) - ∂y(uv) +  forcing
        v_rhs = - ∂x(uv) - ∂y(vv)
        p = get_pressure(u_rhs, v_rhs, ∂x, ∂y, Δ⁻¹)
        u_rhs -= ∂x(p)
        v_rhs -= ∂y(p)
        return [u_rhs, v_rhs]
    end
end

rhs = closure_rhs(∂x, ∂y, Δ⁻¹, P, iP, forcing)

function closure_lhs(Δ, Re, Δt)
    function lhs()
        u_lhs = 1 + -(Δt / Re) * Δ
        return [u_lhs, u_lhs]
    end
end

lhs = closure_lhs(Δ, Re, Δt)
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

##
contourf(x[:], y[:], real.(forcing)', linewidth = 0, color = :thermometer)
##
u = @. 1 * sin(x) * sin(y) + 0im
v = @. 1 * cos(x) * cos(y) + 0im
û = copy(u)
v̂ = copy(v)
mul!(û, P, u)
mul!(v̂, P, v)
# û, v̂ = copy(good_start)
state = [û, v̂]
## Checking Incompressibility
rhs_tmp = rhs(state)
maximum(abs.(∂x(rhs_tmp[1]) + ∂y(rhs_tmp[2])))
##
for i in 1:20000
    imex_step!(state, lhs_operator, rhs, Δt)
    if (i%1000) == 0
        println(i)
        ω̂ = ∂x(v̂) - ∂y(û)
        mul!(u, iP, ω̂)
        local p1 = contourf(x[:],y[:], real.(u)', linewidth = 0, color = :thermometer)
        display(p1)
        println("checking conergence to steady state")
        println(norm(rhs(state) + (1/Re) * Δ.(state)) / (norm(rhs(state))))
        sleep(0.01)
    end 
end
# good_start = copy(state)
## check steady state
# include(pwd() * "/sandbox/two_d_navier_stokes.jl")