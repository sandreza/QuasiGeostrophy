using QuasiGeostrophy, Plots, FFTW, BenchmarkTools
using LinearAlgebra
import Plots: plot
import QuasiGeostrophy: compute


##

Ωxy = Torus(0,2π) × Torus(0,2π)
Nx = 2^8; Ny = 2^8;
fourier_grid = create_grid((Nx, Ny), Ωxy)
x, y = fourier_grid.grid
fourier_transform = Transform(fourier_grid)

fmd  = FourierMetaData("ϕ ", fourier_grid, fourier_transform)
fmd1 = FourierMetaData("ϕ1", fourier_grid, fourier_transform)
fmd2 = FourierMetaData("ϕ2", fourier_grid, fourier_transform)
fmd3 = FourierMetaData("ϕ3", fourier_grid, fourier_transform)
f1 = @. sin(x) + 0im * y
f2 = @. sin(y) + 0im * x
f3 = @. sin(x) * sin(y)
f1 = fourier_transform.forward * f1
f2 = fourier_transform.forward * f2
f3 = fourier_transform.forward * f3
ϕ  = FourierField(f1, fmd)
ϕ1 = FourierField(f1, fmd1)
ϕ2 = FourierField(f2, fmd2)
ϕ3 = FourierField(f3, fmd3)

##
@btime forward(ϕ);
fwd = fourier_transform.forward
dd = ϕ.data
@btime fwd * dd;

##
# Define Closed Operations for FourierField

for unary_operator in unary_operators
    b_symbol = Meta.parse.(unary_operator[2]) #broadcast
    @eval import Base: $b_symbol
    @eval function $b_symbol(field1::FourierField)
        data = broadcast($b_symbol, field1.data)
        metadata  = field1.metadata
        symbname = string($b_symbol)
        name = symbname * "(" * field1.metadata.name * ")"
        fmd = FourierMetaData(name, metadata.grid, metadata.transform)
        FourierField(data, fmd )
    end
end
##

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
∂x = FourierDerivative(im .* kx)
∂y = FourierDerivative(im .* ky)

function (∇::FourierDerivative)(ϕ::FourierField)
    return FourierField(∇.op .* ϕ.data, ϕ.metadata)
end
## Perhaps Define Operator Object

##
dmd = DerivativeMetaData(FourierDerivative(im .* kx), "x")
∂x = Operator(nothing, dmd)
∂x(f_ϕ1)
dmd = DerivativeMetaData(FourierDerivative(im .* ky), "y")
∂y = Operator(nothing, dmd)
∂ʸϕ1 = 2 * ∂y(f_ϕ2)
evaluate(∂ʸϕ1)

plot(compute(∂ʸϕ1))


a1 = compute(∂x)^2
a2 = compute(∂x^2)
norm(a1.op - a2.op)

Δ1 = compute(∂x)^2 + compute(∂y)^2
Δ2 = compute(∂x^2+∂y^2)
norm(Δ1.op - Δ2.op)





