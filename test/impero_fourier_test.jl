using QuasiGeostrophy, FFTW, Test
using LinearAlgebra
include(pwd() * "/test/test_utils.jl")
const tol = 1e1

## Define 1D Test
Ω = S¹(0, 2π) 
Nx = 2^4; 
fourier_grid = FourierGrid(Nx, Ω)
fieldnames = ("ϕ1", "ϕ2", "ϕ3", "ϕ4", "ϕ5")
ϕ1, ϕ2, ϕ3, ϕ4, ϕ5 = create_fields(names = fieldnames, grid = fourier_grid)
∂x = create_operators(fourier_grid)
# initialize FourierField data
x = fourier_grid.grid[1]
ϕ1(sin.(x)); ϕ2(cos.(x));
ϕ3( @. sin(x) * cos(x) );
ϕ4( @. sin(x) + cos(x) );
ϕ5( @. cos(x)^2 - sin(x)^2 );
# Wrap Around Impero Object
@wrapper ϕ1=ϕ1 ϕ2=ϕ2 ϕ3=ϕ3 ϕ4=ϕ4 ϕ5=ϕ5
Δ  = Operator(∂x^2)
∂x = Operator(∂x)

@testset "1D Algebra Test" begin
    println("checking ϕ1 * ϕ2 - ϕ3")
    @test norm(compute(ϕ1 * ϕ2 - ϕ3))/norm(compute(ϕ3)) < eps(tol)
    println("checking ϕ1 + ϕ2 - ϕ4")
    @test norm(compute(ϕ1 + ϕ2 - ϕ4))/norm(compute(ϕ4)) < eps(tol)
    println("checking ϕ1 + ϕ1 - 2 * ϕ1")
    @test norm(compute(ϕ1 + ϕ1- 2*ϕ1))/norm(compute(ϕ1)) < eps(tol)
end

@testset "1D Calculus Test" begin
    println("checking ∂x(ϕ1) - ϕ2")
    @test norm(compute(∂x(ϕ1) - ϕ2))/norm(compute(ϕ2)) < eps(tol)
    println("checking Δ(ϕ1) + ϕ1")
    @test norm(compute(Δ(ϕ1) + ϕ1))/norm(compute(ϕ1)) < eps(tol*10)
end
