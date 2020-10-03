export Operator

struct Operator{𝒮, 𝒯} <: AbstractExpression
    operand::𝒮
    metadata::𝒯
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

struct DerivativeMetaData{𝒪, 𝒟}
    operation::𝒪
    direction::𝒟
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