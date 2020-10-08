using QuasiGeostrophy, Plots, FFTW, BenchmarkTools
using LinearAlgebra
import Plots: plot
import QuasiGeostrophy: compute


##
Ωxy = S¹(0,2π) × S¹(0,2π)
Nx = 2^8; Ny = 2^8;
fourier_grid = create_grid((Nx, Ny), Ωxy)
x, y = fourier_grid.grid
fourier_transform = Transform(fourier_grid)

names = ("ϕ1", "ϕ2", "ϕ3", "ϕ4")
create_fields(names = names, 
              grid = fourier_grid,
              transform = fourier_transform)
f1 = @. sin(x) + 0im * y
f2 = @. sin(y) + 0im * x
f3 = @. sin(x) * sin(y) # product
f4 = @. sin(x) + sin(y) # sum
f1 = fourier_transform.forward * f1
f2 = fourier_transform.forward * f2
f3 = fourier_transform.forward * f3
f4 = fourier_transform.forward * f4
ϕ1.data .= f1
ϕ2.data .= f2
ϕ3.data .= f3
ϕ4.data .= f4
## Check Algebra
norm((ϕ1 * ϕ2 - ϕ3).data)/norm((ϕ3).data)

f_ϕ1 = Field(ϕ1, BasicMetaData("ϕ1"))
f_ϕ2 = Field(ϕ2, BasicMetaData("ϕ2"))
tt =  2 * f_ϕ1 + f_ϕ2 * f_ϕ1 + 2 + tanh(f_ϕ1) + 2*2
compute(a::FourierField) = a
compute(tt)
evaluate(tt)

## Check Calculus
kx, ky = fourier_grid.wavenumbers
∂x = FourierOperator(im .* kx, FourierOperatorMetaData("∂x"))
∂y = FourierOperator(im .* ky, FourierOperatorMetaData("∂y"))

∂x(ϕ)

∂x(ϕ) + 1 + (∂x^2 + ∂y^2)(ϕ)
##
dmd = DerivativeMetaData(FourierOperator(im .* kx, FourierOperatorMetaData("∂x")), "x")
∂x = Operator(nothing, dmd)
∂x(f_ϕ1)
dmd = DerivativeMetaData(FourierOperator(im .* ky, FourierOperatorMetaData("∂y")), "y")
∂y = Operator(nothing, dmd)
∂ʸϕ1 = ∂y(f_ϕ2) * ∂x(f_ϕ1)
evaluate(∂ʸϕ1)
tmp = compute(∂ʸϕ1)
p1 = plot(tmp)


a1 = compute(∂x)^2
a2 = compute(∂x^2)
norm(a1.op - a2.op)

Δ1 = compute(∂x)^2 + compute(∂y)^2
Δ2 = compute(∂x^2+∂y^2)
norm(Δ1.op - Δ2.op)
