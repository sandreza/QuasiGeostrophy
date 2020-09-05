using QuasiGeostrophy, LinearAlgebra

N = 8
x = fourier_nodes(N, a = 0, b = 2π)
y = sin.(x)
plot(x,y)
