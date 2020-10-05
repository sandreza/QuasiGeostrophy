using QuasiGeostrophy, FFTW, Test
using LinearAlgebra
include(pwd() * "/test/test_utils.jl")
const tol = 1e1

## Define 1D Test
Ω = Torus(0, 2π) 
Nx = 2^4; 
fourier_grid = create_grid(Nx, Ω)
fieldnames = ("ϕ1", "ϕ2", "ϕ3", "ϕ4", "ϕ5")
create_fields(names = fieldnames, grid = fourier_grid)
# initialize fields with nontrivial data
x = fourier_grid.grid[1]
ϕ1(sin.(x)); ϕ2(cos.(x));
ϕ3( @. sin(x) * cos(x) );
ϕ4( @. sin(x) + cos(x) );
ϕ5( @. cos(x)^2 - sin(x)^2 );
          
@testset "1D Field Algebra Test" begin
    println("checking ϕ1 * ϕ2 - ϕ3")
    @test norm(ϕ1 * ϕ2 - ϕ3)/norm(ϕ3) < eps(tol)
    println("checking ϕ1 + ϕ2 - ϕ4")
    @test norm(ϕ1 + ϕ2 - ϕ4)/norm(ϕ4) < eps(tol)
end
println(" ")
create_operators(fourier_grid)

@testset "1D Field Calculus Test" begin
    println("checking ∂x(ϕ1)  - ϕ2")
    @test norm(∂x(ϕ1)  - ϕ2)/norm(ϕ2) < eps(tol)
    println("checking ∂x(ϕ3)  - ϕ5")
    @test norm(∂x(ϕ3)  - ϕ5)/norm(ϕ5) < eps(tol)
end

## Define 2D Test
Ω = Torus(0,2π) × Torus(0,2π)
Nx = 2^4; Ny = 2^4;
fourier_grid = create_grid((Nx, Ny), Ω)

fieldnames = ("ϕ1", "ϕ2", "ϕ3", "ϕ4", "ϕ5")
create_fields(names = fieldnames, grid = fourier_grid)

x, y = fourier_grid.grid
ϕ1(sin.(x) .+ 0 .* y)
ϕ2(sin.(y) .+ 0 .* x)
ϕ3(sin.(x) .* sin.(y))
ϕ4(sin.(x) .+ sin.(y))
ϕ5(cos.(x) .* sin.(y))

ϕ4 = FourierField(f4, fmd4)

@testset "2D Field Algebra Test" begin
    @test norm(ϕ1 * ϕ2 - ϕ3)/norm(ϕ3) < eps(tol)
    @test norm(ϕ1 + ϕ2 - ϕ4)/norm(ϕ4) < eps(tol)
end

println(" ")
create_operators(fourier_grid)

@testset "2D Field Calculus Test" begin
    println("checking ∂y(ϕ1) ")
    @test norm(∂y(ϕ1)) < eps(tol)
    println("checking ∂x(ϕ3) - ϕ5 ")
    @test norm( ∂x(ϕ3) - ϕ5 )/norm(ϕ5) < eps(tol)
end

