using QuasiGeostrophy, LinearAlgebra, Test, FFTW, BenchmarkTools

# boiler plate definitions
Nx = 2^8
Ny = 2^8
ax = 0; bx = 2π; ay = 0; by = 4π
gx  = fourier_nodes(Nx, a = ax, b = bx)
gy  = fourier_nodes(Ny, a = ay, b = by)
kx  = fourier_wavenumbers(Nx, L = bx - ax)
ky  = fourier_wavenumbers(Ny, L = by - ay)
x = reshape(gx, (Nx, 1))
y = reshape(gy, (1, Ny))
kx = reshape(kx, (Nx, 1))
ky = reshape(ky, (1, Ny))
f = @. 2 * sin(x) + sin(2*y) + 0im
∂ˣf = @. 2 * cos(x) + 0im * y
∂ʸf = @. 2 * cos(2*y) + 0im * x
Δf = @. - 2 * sin(x) - 4 * sin(2*y)
f̂ = copy(f)
FFTW.set_num_threads(Threads.nthreads())
P = plan_fft(f) # typeof(iP) <: AbstractFFTs.ScaledPlan
iP = plan_ifft(f)
mul!(f̂, P, f)
∂x = FourierDerivative(im .* kx)
∂y = FourierDerivative(im .* ky)
Δ = ∂x^2 + ∂y^2
∫𝒢dxdy = inv(Δ)

@testset "FourierTests" begin
    tolerance = eps(1.0)*1e4
    bool = norm(∂ˣf - iP * ∂x(f̂)) ./ norm(∂ˣf ) < tolerance
    @test bool 
    bool = norm(∂ʸf - iP * ∂y(f̂)) ./ norm(∂ʸf ) < tolerance
    @test bool 
    bool = norm(Δf  - iP * Δ(f̂) ) ./ norm(Δf) < tolerance
    @test bool 
    bool = norm(f - iP * ∫𝒢dxdy(Δ(f̂))) ./ norm(f) < tolerance
    @test bool
end
