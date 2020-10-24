export FourierOperator, FourierOperatorMetaData
export filter_convolve, convolve, box_filter

import Base: *, ^, +, inv, - 
# Derivative Struct and Operations
struct FourierOperator{T, S}
    op::T
    metadata::S
end

struct FourierOperatorMetaData{S}
    name::S
end

FourierOperator(a) = FourierOperator(a, nothing)

function *(∂x::FourierOperator, ϕ::AbstractArray)
    return ∂x.op .* ϕ
end

function (p::FourierOperator)(ϕ::AbstractArray)
    return *(p, ϕ)
end

## Operator Algebra
function ^(p::FourierOperator, α::Number)
    return FourierOperator(p.op .^(α))
end

function *(p::FourierOperator, α::Number)
    return FourierOperator(p.op.*(α))
end
*(q::Number, p::FourierOperator) = *(p, q)


function +(p::FourierOperator, q::FourierOperator)
    return FourierOperator(p.op .+ q.op)
end

function *(p::FourierOperator, q::FourierOperator)
    return FourierOperator(p.op .* q.op)
end

function +(p::FourierOperator, q::Number)
    return FourierOperator(p.op .+ q)
end

function -(p::FourierOperator)
    return FourierOperator( -(p.op))
end

-(p::FourierOperator, q::FourierOperator) = p + (-q)
+(q::Number, p::FourierOperator) = +(p, q)

## with names
function ^(p::FourierOperator{𝒮, 𝒯}, α::Number) where
    {𝒮, 𝒯 <: FourierOperatorMetaData}
    name = p.metadata.name * "^" * string(α)
    fomd = FourierOperatorMetaData(name)
    return FourierOperator(p.op .^(α), fomd)
end

function *(p::FourierOperator{𝒮, 𝒯}, α::Number) where
    {𝒮, 𝒯 <: FourierOperatorMetaData}
    name = string(α) * "*" * p.metadata.name
    fomd = FourierOperatorMetaData(name)
    return FourierOperator(p.op .* α, fomd)
end

function +(p::FourierOperator{𝒮, 𝒯}, q::FourierOperator{𝒮, 𝒯}) where
    {𝒮, 𝒯 <: FourierOperatorMetaData}
    name = "(" * p.metadata.name * "+" * q.metadata.name * ")"
    fomd = FourierOperatorMetaData(name)
    return FourierOperator(p.op .+ q.op, fomd)
end

function *(p::FourierOperator{𝒮, 𝒯}, q::FourierOperator{𝒮, 𝒯}) where
    {𝒮, 𝒯 <: FourierOperatorMetaData}
    name = "(" * p.metadata.name * "*" * q.metadata.name * ")"
    fomd = FourierOperatorMetaData(name)
    return FourierOperator(p.op .* q.op, fomd)
end

function +(p::FourierOperator{𝒮, 𝒯}, q::Number) where
    {𝒮, 𝒯 <: FourierOperatorMetaData}
    name = "(" * p.metadata.name * "+" * string(q) * ")"
    fomd = FourierOperatorMetaData(name)
    return FourierOperator(p.op .+ q, fomd)
end

function -(p::FourierOperator{𝒮, 𝒯}) where
    {𝒮, 𝒯 <: FourierOperatorMetaData}
    name = "-(" * p.metadata.name * ")"
    fomd = FourierOperatorMetaData(name)
    return FourierOperator( -(p.op), fomd)
end

"""
function inv(a::FourierOperator)

# Description
- Takes the inverse of a FourierOperator. If the inverse is zero,
the it returns zero by defualt

# Argument
- `L`: A FourierOperator Object

# Return
- `L⁻¹`: A FourierOperator that represents the inverse of the original
"""
function inv(a::FourierOperator)
    inv_op = 1 ./ a.op 
    @inbounds for i in eachindex(inv_op)
        if abs(inv_op[i]) == Inf
            inv_op[i] = 0.0
        elseif isnan(norm(inv_op[i]))
            inv_op[i] = 0.0
        end
    end
    return FourierOperator(inv_op)
end

function inv(a::FourierOperator{𝒮, 𝒯}) where
    {𝒮, 𝒯 <: FourierOperatorMetaData}
    name = a.metadata.name * "⁻¹"
    fomd = FourierOperatorMetaData(name)
    inv_op = 1 ./ a.op 
    @inbounds for i in eachindex(inv_op)
        if abs(inv_op[i]) == Inf
            inv_op[i] = 0.0
        elseif isnan(norm(inv_op[i]))
            inv_op[i] = 0.0
        end
    end
    return FourierOperator(inv_op, fomd)
end

# Filters
function box_filter(û; pm = 6)
    n = length(û)
    mid = div(n, 2) + 1
    pm  = div(n, pm)
    u = copy(û)
    u[mid - pm : mid + pm] .= 0.0
    return u
end

# Convolution operations
function convolve(û, v̂, P, iP)
    u = iP * û
    v = iP * v̂
    return P * (u .* v)
end

function filter_convolve(û, v̂, P, iP; pm = 6)
    u = box_filter(û, pm = pm)
    v = box_filter(v̂, pm = pm)
    u = iP * u
    v = iP * v
    return P * (u .* v)
end

function convolve!(u, v, w, ŵ, û, v̂, P, iP)
    mul!(u , iP , û)
    mul!(v , iP , v̂)
    @inbounds for i in eachindex(w)
        w[i] = u[i] * v[i]
    end
    mul!(ŵ, P, w)
    return nothing
end

## 
function Base.show(io::IO, O::FourierOperator{𝒮, 𝒯}) where 
    {𝒮, 𝒯 <: FourierOperatorMetaData}
    printstyled(io, O.metadata.name, color = 159)
end