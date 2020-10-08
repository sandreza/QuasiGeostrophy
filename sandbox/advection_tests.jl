using Plots, QuasiGeostrophy
# 1D Tests
Ω = S¹(0,2π)
Nx = 2^8; 
grid = FourierGrid(Nx, Ω)
ϕ, ϕ2 = create_fields(names = ("ϕ", "ϕ2"), grid = grid)
∂x = create_operators(grid)
#  initialize fields with nontrivial data
#  Heuns Method
x = grid.grid[1]

scale = round(Int, 1/(x[2]-x[1]))
ν = 10000
κ = 0e-4
iv = @. tanh( ν*(x - π/2)) * tanh(ν*(x - 3π/2))
ϕ(iv)
ϕ2(iv)
plot(ϕ)
Δt = 0.01 / scale * 2π
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

Ω = S¹(0,2π) × S¹(0,2π)
Nx = Ny =  2^6; 
grid = FourierGrid((Nx, Ny), Ω)
fieldnames = ("ϕ", "ϕ2", "u", "v", "ψ")
ϕ, ϕ2, u, v, ψ = create_fields(names = fieldnames, grid = grid)
∂x, ∂y =  create_operators(grid)
# initialize fields with nontrivial data
#  Heuns Method
x, y = grid.grid

scale = round(Int, 1/(x[2]-x[1]))
ν = 10000
κ = 0e-2
iv = @. tanh( ν*(x - π/2)) * tanh(ν*(x - 3π/2) ) * tanh( ν*(y - π/2)) * tanh(ν*(y - 3π/2) )
# iv = Reshape(iv, (Nx,Ny))
ϕ(iv)
ϕ2(iv)
plot(ϕ)
Δt = 0.01 / scale * 2π
t = randn(1) .* 0
u( (x .+ y) .* 0 .+1 )
v( (x .+ y) .* 0 .+1 )
ψ(@. sin(x/2) * sin(y/2))
# u = ∂y(ψ)
# v = - ∂x(ψ)
norm(∂x(u) + ∂y(v))
for i in 1:100*scale
    Δ = (∂x^2 + ∂y^2)^1
    f = u * ∂x(ϕ) + v * ∂y(ϕ) + κ * Δ(ϕ)
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
    f =  -u * ∂x(ϕ) + -v * ∂y(  ϕ) + κ * Δ(ϕ)
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