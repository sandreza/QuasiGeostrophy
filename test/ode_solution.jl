using QuasiGeostrophy, LinearAlgebra, Test
# Polynomial Order
n = 16

# Domain [a b]
Ω = [-1.0 1.0]
b = Ω[end]
a = Ω[1]
L = b - a

D, cx = cheb(n)

x =  (b-a) .* (cx .+ 1) ./ 2 .+ a

∂ˣ = 2 / L .* D

# Equation Dy = y, y[x=-1] = 1
A = ∂ˣ - I

# Boundary Conditions
A[end, :] .= 0.0
A[end, end] = 1.0

# Solve
rhs = 0.0 * x
rhs[end] = 1.0
y = A \ rhs

# Exact solution
solution = exp.( x .+ 1)

@testset "ODE Test 1" begin
    bool = norm(y - solution) < eps(1.0) * 100
    @test bool
end


# Polynomial Order
n = 32

# Domain [a b]
Ω = [-1.0 0.0]
b = Ω[end]
a = Ω[1]
L = b - a

D, cx = cheb(n)

x =  (b-a) .* (cx .+ 1) ./ 2 .+ a

∂ˣ = 2 / L .* D

# Equation Dy = y, y[x=-1] = 1
A = ∂ˣ - I

# Boundary Conditions
A[end, :] .= 0.0
A[end, end] = 1.0

# Solve
rhs = 0.0 * x
rhs[end] = 1.0
y = A \ rhs

# Exact solution
solution = exp.( x .+ 1)

@testset "ODE Test 2" begin
    bool = norm(y - solution) < eps(1.0) * 100
    @test bool
end
