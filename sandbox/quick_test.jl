using QuasiGeostrophy, LinearAlgebra, Plots, BenchmarkTools
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
luA = lu(A)
y = A \ rhs

# Exact solution
solution = exp.( x .+ 1)

plot(x, solution, label = "analytic", xlabel = "x", ylabel = "solution")
scatter!(x, y, label = "numerical")
plot!(legend = :topleft)

##
# performance check
n = m = 32
A = lu(randn(m*n,m*n))
b = ones(m*n)
x = copy(b)
@btime x = A \ b
@btime ldiv!(x, A, b)
@btime ldiv!(A, b)
