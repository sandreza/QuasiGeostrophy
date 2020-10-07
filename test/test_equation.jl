using QuasiGeostrophy, Test
include(pwd() * "/test/test_utils.jl")
@wrapper u=1 σ=1

@macroexpand  equ = @to_equation σ=u 
equ = @to_equation σ=u
@pde_system pde_system = [
    σ = ∂x(u),
    ∂t(u) = -∂x(u * u - ∂x(σ)),
]

@testset "Basic Equation Test" begin
    @pde_system pde_system = [
        σ = ∂x(u),
        ∂t(u)= -∂x(u * u - ∂x(σ)),
    ]
    @test pde_system[1].rhs == ∂x(u)
end
