using QuasiGeostrophy, LinearAlgebra, Plots, FFTW, BenchmarkTools

N = 8^3
x = fourier_nodes(N, a = 0, b = 2π)
y = sin.(x)
plot(x,y)
tmp = y .+ 0im
tmp2 = y .+ 0im

FFTW.set_num_threads(1)
P = plan_fft(x)
@btime mul!(tmp, P, tmp2)

##
N = 2^8
g  = fourier_nodes(N, a = 0, b = 2π)
x = reshape(g, (N, 1))
y = reshape(g, (1, N))
f = @. sin(x) + sin(y) + 0im
f̂ = copy(f)
FFTW.set_num_threads(6)
P = plan_fft(f)
@btime mul!(f̂, P, f)
##

