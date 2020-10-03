using QuasiGeostrophy, Plots, FFTW, BenchmarkTools
import Plots: plot
import QuasiGeostrophy: compute
import Base: * 
include(pwd() * "/test/test_utils.jl")

struct FourierField{D,S}
    data::D
    metadata::S
end

struct FourierMetaData{𝒩, 𝒢, 𝒯} <: AbstractMetaData 
    name::𝒩
    grid::𝒢
    transform::𝒯
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
        plot(ϕ.metadata.grid.grid[1], ϕ.data)
        plot!(xlabel = "x")
        plot!(ylabel = ϕ.metadata.name)
    elseif dims == 2
        x = ϕ.metadata.grid.grid[1][:]
        y = ϕ.metadata.grid.grid[2][:]
        contourf(x, y, ϕ.data')
        contourf!(xlabel = "x")
        contourf!(ylabel = "y")
    else
        print("Plotting is not supported for fields ")
        print("with dimensions greater ≥ 3")
    end
end
##

Ωxy = Torus(0,2π) × Torus(0,2π)
Nx = 2^8; Ny = 2^8;
fourier_grid = create_grid((Nx, Ny), Ωxy)
x, y = fourier_grid.grid
kx, ky = fourier_grid.wavenumbers
fourier_transform = Transform(fourier_grid)

fmd = FourierMetaData("ϕ", fourier_grid, fourier_transform)
f1 = @. sin(x) + 0im * y
f2 = @. sin(y) + 0im * x
f3 = @. sin(x) * sin(y)
f1 = fourier_transform.forward * f1
f2 = fourier_transform.forward * f2
f3 = fourier_transform.forward * f3
ϕ = FourierField(f1, fmd)
ϕ1 = FourierField(f1, fmd)
ϕ2 = FourierField(f2, fmd)
ϕ3 = FourierField(f3, fmd)

##
@btime forward(ϕ);
fwd = fourier_transform.forward
dd = ϕ.data
@btime fwd * dd;

##
for unary_operator in unary_operators
    b_symbol = Meta.parse.(unary_operator[2]) #broadcast
    @eval import Base: $b_symbol
    @eval $b_symbol(field1::FourierField) where {𝒯} = FourierField(broadcast($b_symbol, field1.data), field1.metadata)
end

for binary_operator in [binary_operators..., ["Negative", "-"]]
    b_symbol = Meta.parse.(binary_operator[2]) #broadcast
    @eval import Base: $b_symbol
    @eval $b_symbol(field1::FourierField, field2::FourierField) = FourierField(broadcast($b_symbol, field1.data, field2.data), field1.metadata)
    @eval $b_symbol(field1::FourierField, field2::𝒮) where {𝒮} =  FourierField(broadcast($b_symbol,field1.data, field2), field1.metadata)
    @eval $b_symbol(field1::𝒯, field2::FourierField) where {𝒯} =  FourierField(broadcast($b_symbol, field1, field2.data), field2.metadata)
    # otherwise there is a method error, data wrapper makes it a closed system
    #@eval $b_symbol(field1::AbstractData, field2::𝒮) where {𝒮 <: Number} = Data(broadcast($b_symbol,field1.data, field2))
    #@eval $b_symbol(field1::𝒯, field2::AbstractData) where {𝒯 <: Number} = Data(broadcast($b_symbol, field1, field2.data))
end
# overwrite multiplication

function *(f̂::FourierField, ĝ::FourierField)
    fwd = f̂.metadata.transform.forward
    bwd = f̂.metadata.transform.forward
    f = bwd * f̂.data
    g = bwd * ĝ.data
    fg = broadcast(*, f, g)
    return FourierField(fwd * fg, f̂.metadata)
end
##
ϕ1 * ϕ2 - ϕ3