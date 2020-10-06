using Plots, QuasiGeostrophy, Test
include(pwd() * "/test/test_utils.jl")
# 1D Tests
Ω = Torus(0,2π) 
Nx = 2^8; 
fourier_grid = create_grid(Nx, Ω)
fieldnames = ("ϕ", "ϕ2")
create_fields(names = fieldnames, grid = fourier_grid)
create_operators(fourier_grid)
# initialize fields with nontrivial data
#  Heuns Method
x = fourier_grid.grid[1]

scale = round(Int, 1/(x[2]-x[1]))
ν = 10000
κ = 0e-4
iv = @. tanh( ν*(x - π/2)) * tanh(ν*(x - 3π/2))
ϕ(iv)
ϕ2(iv)
plot(ϕ)
Δt = 0.01 / scale * 2π
push!(Δt_list_2, Δt)
t = randn(1) .* 0
for i in 1:100*scale
    Δ = ∂x^2
    f = ∂x(ϕ) + κ * Δ(ϕ)
    nϕ = ϕ + Δt * f
    nϕ = ϕ + Δt/2 * ( f +  ∂x(nϕ)  + κ * Δ(nϕ))
    ϕ.data .= nϕ.data
    if i%(scale) == 0
        display(plot(ϕ))
        sleep(0.01)
    end
    t[1] += Δt
end
norm(ϕ-ϕ2)/norm(ϕ2)

##
# 2D Test

Ω = Torus(0,2π) × Torus(0,2π)
Nx = Ny =  2^6; 
fourier_grid = create_grid((Nx, Ny), Ω)
fieldnames = ("ϕ", "ϕ2", "u", "v", "ψ")
create_fields(names = fieldnames, grid = fourier_grid)
create_operators(fourier_grid)
# initialize fields with nontrivial data
#  Heuns Method
x, y = fourier_grid.grid

scale = round(Int, 1/(x[2]-x[1]))
ν = 10000
κ = 0e-2
iv = @. tanh( ν*(x - π/2)) * tanh(ν*(x - 3π/2) ) * tanh( ν*(y - π/2)) * tanh(ν*(y - 3π/2) )
# iv = Reshape(iv, (Nx,Ny))
ϕ(iv)
ϕ2(iv)
plot(ϕ)
Δt = 0.01 / scale * 2π
push!(Δt_list_2, Δt)
t = randn(1) .* 0
#u( (x .+ y) .* 0 .+1 )
#v( (x .+ y) .* 0 .+1 )
ψ(@. sin(x/2) * sin(y/2))
u = ∂y(ψ)
v = - ∂x(ψ)
norm(∂x(u) + ∂y(v))
for i in 1:100*scale
    Δ = (∂x^2 + ∂y^2)^1
    f =  ∂x( u * ϕ) + ∂y( v *  ϕ) + κ * Δ(ϕ)
    nϕ = ϕ + Δt * f
    nϕ = ϕ + Δt/2 * ( f +  u * ∂x(nϕ) + v * ∂y(nϕ)  + κ * Δ(nϕ))
    ϕ.data .= nϕ.data
    if i%(scale) == 0
        display(plot(ϕ))
        sleep(0.01)
    end
    t[1] += Δt
end
for i in 1:100*scale
    Δ = (∂x^2 + ∂y^2)^1
    f =  ∂x( -u * ϕ) + ∂y( -v *  ϕ) + κ * Δ(ϕ)
    nϕ = ϕ + Δt * f
    nϕ = ϕ + Δt/2 * ( f +  -u * ∂x(nϕ) + -v * ∂y(nϕ)  + κ * Δ(nϕ))
    ϕ.data .= nϕ.data
    if i%(scale) == 0
        display(plot(ϕ))
        sleep(0.01)
    end
    t[1] += Δt
end
norm(ϕ-ϕ2) / norm(ϕ2)