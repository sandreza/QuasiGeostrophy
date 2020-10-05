using Plots, QuasiGeostrophy
include(pwd() * "/test/test_utils.jl")
Ω = Torus(0,2π) 
Nx = 2^4; 
fourier_grid = create_grid(Nx, Ω)
fieldnames = ("ϕ1", "ϕ2")
create_fields(names = fieldnames, grid = fourier_grid)
create_operators(fourier_grid)
# initialize fields with nontrivial data
x = fourier_grid.grid[1]
ϕ1(sin.(x)); ϕ2(sin.(x));
##
scale = 100
const Δt = 0.01 / scale * 2π
t = randn(1) .* 0
for i in 1:100*scale
    nϕ = ϕ1 + Δt * ( ∂x(ϕ1) )
    ϕ1.data .= nϕ.data
    if i%scale == 0
    display(plot(ϕ1))
    sleep(0.01)
    end
    t[1] += Δt
end

norm(ϕ1 - ϕ2) / Δt # first order convergence