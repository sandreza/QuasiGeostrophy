using QuasiGeostrophy, LinearAlgebra, Test, FFTW, BenchmarkTools, Plots

a = 0; b = 22.0
const Δt = 1e-1 / 2
Ω = Torus(a,b) # Domain
Nx = 2^6;
fourier_grid = create_grid(Nx, Ω) # Grid
x = fourier_grid.grid
kx = fourier_grid.wavenumbers
f = @. sin( 4 * 2π * x / (b-a)) + 0im
f̂ = copy(f)
FFTW.set_num_threads(Threads.nthreads())
P = plan_fft(f) # typeof(iP) <: AbstractFFTs.ScaledPlan
iP = plan_ifft(f)
mul!(f̂, P, f)
∂x = FourierDerivative(im .* kx)
argmax(real.((-∂x^2 - ∂x^4).op))

function closure_rhs(∂x, P, iP)
    function closed_rhs(state)
        f̂ = state[1]
        ff = convolve(f̂, f̂, P, iP)
        return [ ∂x( -0.5 .* ff - ∂x(f̂))]
    end
end

rhs = closure_rhs(∂x, P, iP)

function closure_lhs(∂x, Δt)
    function lhs()
        q1_lhs = 1 + Δt * ∂x^4
        return [q1_lhs]
    end
end

lhs = closure_lhs(∂x, Δt)
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
f = @. sin( 8π * x / (b-a)) + 0im
f̂ = copy(f)
mul!(f̂, P, f)
state = [f̂]
nsteps = 40*200
total_array = zeros(length(f), nsteps) 
##
gr(size = (400,400))
for i in 1:nsteps
    imex_step!(state, lhs_operator, rhs, Δt)
    u = state[1]
    total_array[:, i] .= real.(iP * u)
    if (i%100) == 0
        mul!(f, iP, f̂)
        p1 = plot(x, real.(f), ylims = (-2.0,2.0), label = "real space", xlabel = "position", ylabel  = "amplitude")
        p2 = scatter(log.(abs.(f̂))[1:div(length(f̂),2)+1], label = "spectrum", ylims = (-7, 3), xlabel = "wavenumber", ylabel = "spectral amplitude")
        display(plot(p1,p2))
        sleep(0.1)
    end
end
##
contourf(total_array, linewidth = 0, color = :thermometer)