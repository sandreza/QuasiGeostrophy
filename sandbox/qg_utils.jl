import QuasiGeostrophy: compute
struct Wrapper{T, S} <: AbstractExpression
    data::T
    meta_data::S
end
struct WrapperMetaData{T}
    io_name::T
end

function Base.show(io::IO, field::Wrapper{T, S}) where {T <: Char, S}
    color = 230
    printstyled(io, field.data, color = color)
end
function Base.show(io::IO, field::Wrapper{T, S}) where {T, S <: WrapperMetaData}
    color = 230
    printstyled(io, field.meta_data.io_name, color = color)
end

compute(a::Wrapper) = a.data

macro wrapper(expr)
    rewritten_expr = _wrapper(expr)
    return rewritten_expr
end

function _wrapper(expr::Expr)
    symb = expr.args[1]
    val  = expr.args[2]
    if expr.head != :(=)
        println( "@wrapper macro not in proper form")
        println( "must be ")
        println( "@wrapper a=1 b=2 c=3")
        return error()
    end
    string_symb = String(symb)
    new_expr = :($(esc(symb)) =  Wrapper($val, WrapperMetaData($string_symb)))
    return new_expr
end

macro wrapper(exprs...)
    rewritten_exprs = [_wrapper(expr) for expr in exprs]
    return Expr(:block, rewritten_exprs...)
end

function rhs(qg_system)
    qg_system[1].lhs.data.data .= evaluate(qg_system[1].rhs)
    qg_system[2].lhs.data.data .= evaluate(qg_system[2].rhs)
    rhs¹ = evaluate(qg_system[3].rhs)
    rhs² = evaluate(qg_system[4].rhs)
    return [rhs¹, rhs²]
end
function evolve_system(qg_system, Δt; filter = 1)
    data = [qg_system[3].lhs.operand.data.data, qg_system[4].lhs.operand.data.data]
    olddata = copy(data)
    stage1 = rhs(qg_system)
    data .+= Δt .* stage1
    stage2 = rhs(qg_system)
    data .= olddata .+ Δt/2 .* (stage1 + stage2)
    qg_system[3].lhs.operand.data.data .= filter * data[1]
    qg_system[4].lhs.operand.data.data .= filter * data[2]
    return nothing
end
#
function cfl(ψ¹, ψ², ∂x, ∂y)
    F = ψ¹.data.metadata.transform.backward
    x,y = ψ¹.data.metadata.grid.grid
    Δx = x[2] - x[1]
    Δy = y[2] - y[1]
    u1 = F * evaluate(∂y(ψ¹))
    v1 = F * evaluate(∂x(ψ¹))
    u2 = F * evaluate(∂y(ψ²))
    v2 = F * evaluate(∂x(ψ²))
    u1 = maximum(abs.(u1))
    v1 = maximum(abs.(v1))
    u2 = maximum(abs.(u2))
    v2 = maximum(abs.(v2))
    return minimum([Δx / u1, Δx / u2, Δy / v1, Δy / v2]) # gives a Δt
end

# building some cutoff filters
function cutoff(k, cutoff)
    return k > cutoff ? 0.0 : 1.0
end
 