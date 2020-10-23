using QuasiGeostrophy, CUDA, BenchmarkTools, FFTW
# 1D Tests file: include("sandbox/gpu_tests.jl")
Ω = S¹(0,2π) 
Nx = 2^8; 
arraytype = CuArray
grid = FourierGrid(Nx, Ω, arraytype = arraytype)
ϕ, ϕ2 = create_fields(names = ("ϕ", "ϕ2"), grid = grid, arraytype = arraytype)
∂x = create_operators(grid, arraytype = arraytype)
# initialize fields with nontrivial data
x = grid.grid[1]
k = grid.wavenumbers[1]

P = ϕ.metadata.transform.forward
c_P = plan_fft(ϕ.data)
scale = round(Int, 1/(x[2]-x[1]))
ν = 10^4
κ = 1e-4
iv = @. tanh( ν * (x - π/2) ) * tanh( ν * (x - 3π/2) )
ϕ(iv)
ϕ2(iv)
Δ =  FourierOperator(∂x.op .* ∂x.op,FourierOperatorMetaData("Δ"))
Δt = 0.01 / scale * 2π
t = randn(1) .* 0

for i in 1:100*scale
    f = ∂x(ϕ) + κ * Δ(ϕ)
    nϕ = ϕ + Δt * f
    g  = ∂x(nϕ)  + κ * Δ(nϕ)
    nϕ = ϕ + Δt/2 * ( f +  g)
    ϕ.data .= nϕ.data
    t[1] += Δt
end
println(norm(ϕ-ϕ2)/norm(ϕ2))

## 2D Tests
# 2D Test
Ω = S¹(0,2π) × S¹(0,2π)
Nx = Ny =  2^8; 
grid = FourierGrid((Nx, Ny), Ω, arraytype = arraytype)
fieldnames = ("ϕ", "ϕ2", "u", "v", "ψ")
ϕ, ϕ2, u, v, ψ = create_fields(names = fieldnames, grid = grid, arraytype = arraytype)
∂x, ∂y =  create_operators(grid, arraytype = arraytype)
Δ =  FourierOperator(∂x.op .* ∂x.op .+ ∂y.op .* ∂y.op,FourierOperatorMetaData("Δ"))
# initialize fields with nontrivial data
#  Heuns Method
x, y = grid.grid
scale = round(Int, 1/(x[2]-x[1]))
ν = 10000
κ = 1e-4
iv = @. tanh( ν*(x - π/2)) * tanh(ν*(x - 3π/2) ) * tanh( ν*(y - π/2)) * tanh(ν*(y - 3π/2) )
ϕ(iv)
ϕ2(iv)
Δt = 0.01 / scale * 2π
t = randn(1) .* 0
u( (x .+ y) .* 0 .+1 )
v( (x .+ y) .* 0 .+1 )
ψ(@. sin(x/2) * sin(y/2))
norm(∂x(u) + ∂y(v))
for i in 1:100*scale
    f = u * ∂x(ϕ) + v * ∂y(ϕ) + κ * Δ(ϕ)
    nϕ = ϕ + Δt * f
    nϕ = ϕ + Δt/2 * ( f +  u * ∂x(nϕ) + v * ∂y(nϕ)  + κ * Δ(nϕ))
    ϕ.data .= nϕ.data
    t[1] += Δt
end
for i in 1:100*scale
    f =  -u * ∂x(ϕ) + -v * ∂y(ϕ) + κ * Δ(ϕ)
    nϕ = ϕ + Δt * f
    nϕ = ϕ + Δt/2 * ( f +  -u * ∂x(nϕ) + -v * ∂y(nϕ)  + κ * Δ(nϕ))
    ϕ.data .= nϕ.data
    t[1] += Δt
end
println(norm(ϕ-ϕ2) / norm(ϕ2))
