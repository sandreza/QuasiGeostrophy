export FourierField, FourierMetaData, Transform
export plot, forward, backward

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