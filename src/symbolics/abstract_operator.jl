export Operator, DerivativeMetaData, OperatorMetaData
export compute

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

struct OperatorMetaData{𝒪, 𝒩}
    operation::𝒪
    name::𝒩
end

function Base.show(io::IO, o::Operator{S,T}) where
    {S <: Nothing, T <: OperatorMetaData}
    name = o.metadata.name
    printstyled(io, name, color = 14 )
end

function Base.show(io::IO, o::Operator{S,T}) where 
    {S <: AbstractExpression, T <: OperatorMetaData}
    name = o.metadata.name
    printstyled(io, name, "(",  color = 14 )
    print(o.operand)
    printstyled(io, ")",  color = 14 )
end

to_expr(t::Operator) = Expr(:call, t, to_expr(t.operand))