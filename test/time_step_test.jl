using Plots, QuasiGeostrophy, Test
include(pwd() * "/test/test_utils.jl")
Ω = Torus(0,2π) 
Nx = 2^3; 
plot_flag = false
plot_flag2 = false
fourier_grid = create_grid(Nx, Ω)
fieldnames = ("ϕ_numerical", "ϕ_exact")
create_fields(names = fieldnames, grid = fourier_grid)
create_operators(fourier_grid)
# initialize fields with nontrivial data
x = fourier_grid.grid[1]
ϕ_exact(sin.(x));
# Forward Euler
Δt_list = []
error_list = []
for scale in [1, 10^1, 10^2, 10^3, 10^4]
    ϕ_numerical(sin.(x))
    local Δt = 0.1 / scale * 2π
    push!(Δt_list, Δt)
    t = randn(1) .* 0
    for i in 1:10*scale
        f = ∂x(ϕ_numerical)
        nϕ = ϕ_numerical + Δt * f
        ϕ_numerical.data .= nϕ.data
        if i%(scale) == 0
            if plot_flag
                display(plot(ϕ_numerical))
                sleep(0.01)
            end
        end
        t[1] += Δt
    end
    error = norm(ϕ_numerical - ϕ_exact)
    push!(error_list, error)
end
#  Heuns Method
Δt_list_2 = []
error_list_2 = []
for scale in [1, 10^1, 10^2, 10^3, 10^4]
    ϕ_numerical(sin.(x))
    local Δt = 0.1 / scale * 2π
    push!(Δt_list_2, Δt)
    t = randn(1) .* 0
    for i in 1:10*scale
        f = ∂x(ϕ_numerical)
        nϕ = ϕ_numerical + Δt * f
        nϕ = ϕ_numerical + Δt/2 * ( f +  ∂x(nϕ) )
        ϕ_numerical.data .= nϕ.data
        if i%(scale) == 0
            if plot_flag
                display(plot(ϕ_numerical))
                sleep(0.01)
            end
        end
        t[1] += Δt
    end
    error = norm(ϕ_numerical - ϕ_exact)
    push!(error_list_2, error)
end
#
prefactor = error_list ./ Δt_list # first order convergence
prefactor2 = error_list_2 ./ Δt_list_2 .^2 # second order convergence
if plot_flag2
    xc = collect(0:0.1:4.5)
    plot(xc, -1 .* xc .+ log10(prefactor[end]), label = "1st order theoretical", xlims = (0, 4.5), linestyle = :dash, linewidth = 2, color = :orange)
    scatter!(-log10.(Δt_list), log10.(error_list), label = "1st order numerical", color = :orange)
    p1 = plot!(xlabel = "-log10(Δt)", ylabel = "log10(error)")
    plot!(xc, -2 .* xc .+ log10(prefactor2[end]), label = "2nd order theoretical", linestyle = :dash, linewidth = 2, color = :green)
    scatter!(-log10.(Δt_list_2), log10.(error_list_2), label = "2nd order numerical", color = :green)
    p2 = plot!(xlabel = "-log10(Δt)", ylabel = "log10(error)")
    display(p2)
end

@testset "TimeStepping Convergence Test" begin
    # forward euler prefactor for this test is about 40
    @test prod(prefactor[2:end] .< 40)
    # Heun's method prefactor for this test is about 12
    @test prod(prefactor2[2:end] .< 12)
end
