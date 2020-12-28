export FourierField, FourierMetaData, Transform
export plot, forward, backward, norm
import LinearAlgebra: norm
struct FourierField{D,S}
    data::D
    metadata::S
end

struct FourierMetaData{𝒩, 𝒢, 𝒯} 
    name::𝒩
    grid::𝒢
    transform::𝒯
end

function Base.show(io::IO, ϕ::FourierField{S,T}) where {S, T <: FourierMetaData}
    printstyled(io, ϕ.metadata.name, color = 128 )
end

function forward(ϕ::FourierField{S,T}) where {S, T <: FourierMetaData}
    ϕ̂ = ϕ.metadata.transform.forward * ϕ.data
    fmd = FourierMetaData(nothing, ϕ.metadata.grid, ϕ.metadata.transform)
    return FourierField(ϕ̂, fmd)
end

function backward(ϕ::FourierField{S,T}) where{S, T <: FourierMetaData}
    ϕ̂ = ϕ.metadata.transform.backward * ϕ.data
    fmd = FourierMetaData(nothing, ϕ.metadata.grid, ϕ.metadata.transform)
    return FourierField(ϕ̂, fmd)
end

# convenience function to initialize array
function (ϕ::FourierField)(a::AbstractArray)
    f1 = ϕ.metadata.transform.forward * (a .+ 0im * (ϕ.data .* false))
    ϕ.data .= f1
    return nothing
end

norm(ϕ::FourierField) = norm(ϕ.data)