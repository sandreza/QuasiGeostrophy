using QuasiGeostrophy
include(pwd() * "/test/test_utils.jl")
import QuasiGeostrophy: @pde_system, @to_equation, _to_equation
@wrapper u=1 σ=1

@macroexpand  equ = @to_equation σ=u 
equ = @to_equation σ=u
@pde_system pde_system = [
    σ = ∂x(u),
    ∂t(u)= -∂x(u * u - ∂x(σ)),
]
