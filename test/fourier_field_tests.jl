using QuasiGeostrophy, FFTW, Test
using LinearAlgebra
include(pwd() * "/test/test_utils.jl")
const tol = 1e1

## Define 1D Test
Ωxy = Torus(0,2π) 
Nx = 2^4; 
fourier_grid = create_grid(Nx, Ωxy)
x = fourier_grid.grid[1]
fourier_transform = Transform(fourier_grid)

names = ("ϕ1", "ϕ2", "ϕ3", "ϕ4")
create_fields(names = names, 
              grid = fourier_grid,
              transform = fourier_transform)
ϕ1(sin.(x)); ϕ2(cos.(x));
ϕ3( @. sin(x) * cos(x) );
ϕ4( @. sin(x) + cos(x) );

@testset "1D Field Test" begin
    @test norm(ϕ1 * ϕ2 - ϕ3)/norm(ϕ3) < eps(tol)
    @test norm(ϕ1 + ϕ2 - ϕ4)/norm(ϕ4) < eps(tol)
end

## Define 2D Test
Ωxy = Torus(0,2π) × Torus(0,2π)
Nx = 2^4; Ny = 2^4;
fourier_grid = create_grid((Nx, Ny), Ωxy)

x, y = fourier_grid.grid
fourier_transform = Transform(fourier_grid)
fmd  = FourierMetaData("ϕ" , fourier_grid, fourier_transform)
fmd1 = FourierMetaData("ϕ1", fourier_grid, fourier_transform)
fmd2 = FourierMetaData("ϕ2", fourier_grid, fourier_transform)
fmd3 = FourierMetaData("ϕ3", fourier_grid, fourier_transform)
fmd4 = FourierMetaData("ϕ4", fourier_grid, fourier_transform)
f1 = @. sin(x) + 0im * y
f2 = @. sin(y) + 0im * x
f3 = @. sin(x) * sin(y) + 0im
f4 = @. sin(x) + sin(y) + 0im
f1 = fourier_transform.forward * f1
f2 = fourier_transform.forward * f2
f3 = fourier_transform.forward * f3
f4 = fourier_transform.forward * f4
ϕ  = FourierField(f1, fmd)
ϕ1 = FourierField(f1, fmd1)
ϕ2 = FourierField(f2, fmd2)
ϕ3 = FourierField(f3, fmd3)
ϕ4 = FourierField(f4, fmd4)

@testset "2D Field Test" begin
    @test norm((ϕ1 * ϕ2 - ϕ3).data)/norm((ϕ3).data) < eps(tol)
    @test norm((ϕ1 + ϕ2 - ϕ4).data)/norm((ϕ4).data) < eps(tol)
end