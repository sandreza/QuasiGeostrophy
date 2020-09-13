export FourierDerivative
export filter_convolve, convolve, box_filter

import Base: *, ^, +, inv, - 
# Derivative Struct and Operations
struct FourierDerivative{T}
    op::T
end

function *(∂x::FourierDerivative, ϕ::AbstractArray)
    return ∂x.op .* ϕ
end

function (p::FourierDerivative)(ϕ::AbstractArray)
    return *(p, ϕ)
end

function ^(p::FourierDerivative, α::Number)
    return FourierDerivative(p.op.^(α))
end

function *(p::FourierDerivative, α::Number)
    return FourierDerivative(p.op.*(α))
end
*(q::Number, p::FourierDerivative) = *(p, q)


function +(p::FourierDerivative, q::FourierDerivative)
    return FourierDerivative(p.op .+ q.op)
end
function +(p::FourierDerivative, q::Number)
    return FourierDerivative(p.op .+ q)
end
function -(p::FourierDerivative)
    return FourierDerivative( -(p.op))
end
-(p::FourierDerivative, q::FourierDerivative) = p + (-q)
+(q::Number, p::FourierDerivative) = +(p, q)

function inv(a::FourierDerivative)
    inv_op = 1 ./ a.op 
    @inbounds for i in eachindex(inv_op)
        if abs(inv_op[i]) == Inf
            inv_op[i] = 0.0
        elseif isnan(norm(inv_op[i]))
            inv_op[i] = 0.0
        end
    end
    return FourierDerivative(inv_op)
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
