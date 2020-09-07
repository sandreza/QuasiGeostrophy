using QuasiGeostrophy, LinearAlgebra, Test, FFTW, BenchmarkTools, Plots

a = 0; b = 2π
Ω = Torus(a,b) # Domain
Nx = 2^8;
fourier_grid = create_grid(Nx, Ω) # Grid
x = fourier_grid.grid
kx = fourier_grid.wavenumbers
f = @. exp(-2 * (b-a) / 3 * (x - (b-a)/2)^2)*(1+0im)
f̂ = copy(f)
FFTW.set_num_threads(Threads.nthreads())
P = plan_fft(f) # typeof(iP) <: AbstractFFTs.ScaledPlan
iP = plan_ifft(f)
mul!(f̂, P, f)
∂x = FourierDerivative(im .* kx)
ν = 0.01

function rhs(f̂, ∂x, P, iP, ν)
    ff = filter_convolve(f̂, f̂, P, iP)
    return -∂x(ff) + ν * ∂x(∂x(f̂))
end
function closure_rhs(∂x, P, iP, ν)
    function closed_rhs(f̂)
        ff = filter_convolve(f̂, f̂, P, iP)
        #ff = convolve(f̂, f̂, P, iP)
        return -∂x(ff) + ν * ∂x(∂x(f̂))
    end
end
old_rhs(f̂) = rhs(f̂, ∂x, P, iP, ν)
new_rhs = closure_rhs(∂x, P, iP, ν)
plot(x, real.(f))
plot(x, real.(iP * ∂x(f̂)))
plot(x, real.(iP * rhs(f̂, ∂x, P, iP, ν)))

function step!(f̂, rhs, Δt)
    rhsf = rhs(f̂)
    @. f̂ += Δt * rhsf
end
function imex_step!(f̂, rhs, Δt)
    ff = filter_convolve(f̂, f̂, P, iP)
    rhsf = -∂x(ff)
    @. f̂ += Δt * rhsf
    lhsf = 1 .- Δt .* ν .* (∂x^2).op
    @. f̂ = f̂ / lhsf
end
##
f = @. exp(-2 * (b-a) / 3 * (x - (b-a)/2)^2)*(1+0im)
f̂ = copy(f)
mul!(f̂, P, f)
for i in 1:2π*200
    local Δt = 0.001
    #step!(f̂, new_rhs, Δt)
    imex_step!(f̂, new_rhs, Δt)
    if (i%100) == 0
        mul!(f, iP, f̂)
        p1 = plot(x, real.(f), ylims = (0,1))
        p2 = scatter(log.(abs.(f̂))[1:div(length(f̂),2)+1], label = "spectrum", ylims = (-7, 3))
        display(plot(p1,p2))
        sleep(0.1)
    end
end
##

plot(x, real.(f), ylims = (0,1))
plot!(x, real.(iP * tmp))
