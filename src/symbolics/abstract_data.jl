export AbstractData, Data, evaluate

abstract type AbstractData end
struct Data{𝒯, 𝒟} <: AbstractData
    data::𝒯
    metadata::𝒟
end

Data(data) = Data(data, nothing)

function Base.show(io::IO, d::Data)
    printstyled(io, d.data, color = 199)
end

for unary_operator in unary_operators
    b_symbol = Meta.parse.(unary_operator[2]) #broadcast
    @eval $b_symbol(field1::AbstractData) where {𝒯} = Data(broadcast($b_symbol, field1.data))
end

for binary_operator in [binary_operators..., ["Negative", "-"]]
    b_symbol = Meta.parse.(binary_operator[2]) #broadcast
    @eval $b_symbol(field1::AbstractData, field2::AbstractData) = Data(broadcast($b_symbol, field1.data, field2.data))
    @eval $b_symbol(field1::AbstractData, field2::𝒮) where {𝒮} =  Data(broadcast($b_symbol,field1.data, field2))
    @eval $b_symbol(field1::𝒯, field2::AbstractData) where {𝒯} = Data(broadcast($b_symbol, field1, field2.data))
    # otherwise there is a method error, data wrapper makes it a closed system
    @eval $b_symbol(field1::AbstractData, field2::𝒮) where {𝒮  <: Number} = Data(broadcast($b_symbol,field1.data, field2))
    @eval $b_symbol(field1::𝒯, field2::AbstractData) where {𝒯 <: Number} = Data(broadcast($b_symbol, field1, field2.data))
end

# Define evaluation functions
compute(d::AbstractData) = d
evaluate(a::AbstractExpression) = compute(a).data
