using QuasiGeostrophy, Plots, FFTW, BenchmarkTools
using LinearAlgebra
import Plots: plot
import QuasiGeostrophy: compute

struct FourierField{D,S}
    data::D
    metadata::S
end

struct FourierMetaData{𝒩, 𝒢, 𝒯} <: AbstractMetaData 
    name::𝒩
    grid::𝒢
    transform::𝒯
end

function Base.show(io::IO, ϕ::FourierField{S,T}) where {S, T <: FourierMetaData}
    printstyled(io, ϕ.metadata.name, color = 128 )
end



struct Transform{ℱ, ℬ}
    forward::ℱ
    backward::ℬ
end

function Transform(𝒢::FourierGrid)
    grid_size = length.(fourier_grid.grid)
    f = randn(grid_size...) .+ 0im
    FFTW.set_num_threads(Threads.nthreads())
    P = plan_fft(f) # typeof(iP) <: AbstractFFTs.ScaledPlan
    iP = plan_ifft(f)
    return Transform(P, iP)
end

function forward(ϕ::FourierField{S,T}) where{S, T<:FourierMetaData}
    ϕ̂ = ϕ.metadata.transform.forward * ϕ.data
    fmd = FourierMetaData(nothing, ϕ.metadata.grid, ϕ.metadata.transform)
    return FourierField(ϕ̂, fmd)
end

function backward(ϕ::FourierField{S,T}) where{S, T<:FourierMetaData}
    ϕ̂ = ϕ.metadata.transform.backward * ϕ.data
    fmd = FourierMetaData(nothing, ϕ.metadata.grid, ϕ.metadata.transform)
    return FourierField(ϕ̂, fmd)
end

function plot(ϕ::FourierField{S, T}) where {S, T <: FourierMetaData}
    dims = length(ϕ.metadata.grid.grid)
    if dims == 1
        x = ϕ.metadata.grid.grid[1][:]
        dd = ϕ.metadata.transform.backward * ϕ.data
        contourf(x, real.(dd))
        plot!(xlabel = "x")
        plot!(ylabel = ϕ.metadata.name)
    elseif dims == 2
        x = ϕ.metadata.grid.grid[1][:]
        y = ϕ.metadata.grid.grid[2][:]
        dd = ϕ.metadata.transform.backward * ϕ.data
        contourf(x, y, real.(dd)')
        contourf!(xlabel = "x")
        contourf!(ylabel = "y")
        contourf!(title =  ϕ.metadata.name)
    else
        print("Plotting is not supported for fields ")
        print("with dimensions greater ≥ 3")
    end
end

plot(ϕ::Field{S, T}) where {S <: FourierField, T} = plot(ϕ.data)
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
for binary_operator in [binary_operators..., ["Negative", "-"]]
    b_symbol = Meta.parse.(binary_operator[2]) #broadcast
    @eval import Base: $b_symbol
    @eval function $b_symbol(field1::FourierField, field2::FourierField)
        data = broadcast($b_symbol, field1.data, field2.data)
        metadata  = field1.metadata
        symbname = string($b_symbol)
        name1 = field1.metadata.name 
        name2 = field2.metadata.name 
        name = "(" * name1 * symbname * name2 * ")"
        fmd = FourierMetaData(name, metadata.grid, metadata.transform)
        return FourierField(data, fmd)
    end
    @eval function $b_symbol(field1::FourierField, field2::𝒮) where {𝒮}
        data = broadcast($b_symbol, field1.data, field2)
        metadata  = field1.metadata
        symbname = string($b_symbol)
        name1 = field1.metadata.name 
        name2 = string(field2)
        name = "(" * name1 * symbname * name2 * ")"
        fmd = FourierMetaData(name, metadata.grid, metadata.transform)
        return FourierField(data, fmd)
    end
    @eval function $b_symbol(field1::𝒯, field2::FourierField) where {𝒯}
        data = broadcast($b_symbol, field1, field2.data)
        metadata  = field2.metadata
        symbname = string($b_symbol)
        name1 = string(field1)
        name2 = field2.metadata.name 
        name = "(" * name1 * symbname * name2 * ")"
        fmd = FourierMetaData(name, metadata.grid, metadata.transform)
        return FourierField(data, fmd)
    end
end
# overwrite multiplication
function *(f̂::FourierField, ĝ::FourierField)
    fwd = f̂.metadata.transform.forward
    bwd = f̂.metadata.transform.backward
    f = bwd * f̂.data
    g = bwd * ĝ.data
    fg = broadcast(*, f, g)
    metadata  = ĝ.metadata
    name1 = f̂.metadata.name
    name2 = ĝ.metadata.name 
    name = "(" * name1 * "*" * name2 * ")"
    fmd = FourierMetaData(name, metadata.grid, metadata.transform)
    return FourierField(fwd * fg, fmd)
end

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

struct Operator{𝒮, 𝒯} <: AbstractExpression
    operand::𝒮
    metadata::𝒯
end
struct DerivativeMetaData{𝒪, 𝒟}
    operation::𝒪
    direction::𝒟
end

function (o::Operator)(expr::AbstractExpression)
    return Operator(expr, o.metadata)
end

function compute(o::Operator)
    return o.metadata.operation(compute(o.operand))
end

function compute(o::Operator{𝒮, 𝒯}) where 
    {𝒮 <: Nothing, 𝒯}
    return compute(o.metadata)
end

function compute(a::DerivativeMetaData{𝒮,𝒯}) where
    {𝒮 <: FourierDerivative, 𝒯}
    return a.operation
end 


function Base.show(io::IO, o::Operator{S,T}) where
    {S <: Nothing, T <: DerivativeMetaData}
    name = Char(0x02202) * o.metadata.direction
    printstyled(io, name, color = 14 )
end

function Base.show(io::IO, o::Operator{S,T}) where 
    {S <: AbstractExpression, T <: DerivativeMetaData}
    name = Char(0x02202) * o.metadata.direction
    printstyled(io, name, "(",  color = 14 )
    print(o.operand)
    printstyled(io, ")",  color = 14 )
end
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





