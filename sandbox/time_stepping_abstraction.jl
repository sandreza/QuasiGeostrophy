using QuasiGeostrophy

abstract type SpatialDiscretization <: AbstractExpression end
abstract type TimeDiscretization <: AbstractExpression end
struct Fourier <: SpatialDiscretization end
struct Time <: TimeDiscretization end

struct AB1 <: TimeDiscretization end
struct BDF1 <: TimeDiscretization end

struct Explicit{T,S} <: TimeDiscretization
    operand::T
    meta_data::S
end

function extract(a, b)
    if typeof(a) <: b
        return a
    else
        return Field(0)
    end
    return nothing
end

##
struct Implicit{T,S} <: TimeDiscretization
    operand::T
    meta_data::S
end

struct AbstractEquation{T,S} <: AbstractExpression
    lhs::T
    rhs::S
end
struct Wrapper{T, S} <: AbstractExpression
    data::T
    meta_data::S
end
# Struct for MetaData
struct MetaData{T}
    io_name::T
end
struct GradientMetaData{T}
    direction::T
end
function Base.show(io::IO, field::Wrapper{T, S}) where {T, S <: MetaData}
    color = 226
    printstyled(io, field.meta_data.io_name, color = color)
end

u = Wrapper("u", MetaData("u"))
σ = Wrapper("σ", MetaData("σ"))
ν = Wrapper(0.001, MetaData("ν"))
∂x(u) = Gradient(u, GradientMetaData("x"))
∂t(u) = Gradient(u, GradientMetaData("t"))
tmp = ∂x(u)
##
function Base.show(io::IO, ∇::Gradient) 
    color = 220
    printstyled(io, "∇(", color = color)
    print(∇.operand)
    printstyled(io, ")", color = color)
end

function Base.show(io::IO, ∇::Gradient) 
    color = 220
    printstyled(io, "∇(", color = color)
    print(∇.operand)
    printstyled(io, ")", color = color)
end

function Base.show(io::IO, ∇::Gradient{T,S}) where {T, S <: GradientMetaData}
    color = 220
    tmp = ∇.metadata.direction
    printstyled(io, "∂", tmp, "(", color = color)
    print(∇.operand)
    printstyled(io, ")", color = color)
end

function Base.show(io::IO, ϕ::Field) 
    color = 220
    print(ϕ.data)
end
##
import Base: ==

function ==(a::AbstractExpression, b::AbstractExpression)
    return AbstractEquation(a,b)
end

function Base.show(io::IO, eq::AbstractEquation) 
    color = 118
    print(io, eq.lhs)
    printstyled(io, "=", color = color)
    print(io, eq.rhs)
end

function Base.show(io::IO, eq::Implicit) 
    print(io, eq.operand)
end
function Base.show(io::IO, eq::Explicit) 
    print(io, eq.operand)
end

function Base.show(io::IO, system::Array{T}) where {T <: AbstractEquation}
    for equation in system
        println(io, equation)
    end
end

∂t(u) == -∂x(u*u) + ν * ∂x(∂x(u)) 

pde_equation = [
    Implicit(∂t(u), BDF1()) == Explicit(-∂x(u*u), AB1()) + Implicit(ν * ∂x(∂x(u)) , BDF1())
]

## Implicit / Explicit grabbing rules.
# Assumes seperation by addition should be okay since quadrature is a 
# linear operation
function get_implicit(a::Implicit)
    return a
end
function get_implicit(a::Explicit)
    return Field(0) # probably need a rule like 0 + field = field
end
function get_implicit(a::Add)
    t1 = get_implicit(a.term1)
    t2 = get_implicit(a.term2)
    return t1 + t2
end
function get_implicit(a::Negative)
    t1 = get_implicit(a.term)
    return t1
end

function get_explicit(a::Explicit)
    return a
end

function get_explicit(a::Implicit)
    return Field(0) # probably need a rule like 0 + field = field
end

function get_explicit(a::Add) 
    t1 = get_explicit(a.term1)
    t2 = get_explicit(a.term2)
    return t1 + t2
end
function get_explicit(a::Negative)
    t1 = get_explicit(a.term)
    return t1
end

function get_implicit(eq::AbstractEquation)
    # place all implicit terms on left-hand side
    implicit = get_implicit(eq.lhs) - get_implicit(eq.rhs)
    return implicit
end
function get_explicit(eq::AbstractEquation)
    # place all explicit terms on left-hand side
    explicit = get_explicit(eq.rhs) - get_explicit(eq.lhs)
    return explicit
end
function get_imex(eq::AbstractEquation)
    implicit = get_implicit(eq)
    explicit = get_explicit(eq)
    return implicit == explicit
end
function get_imex(system::Array{S}) where {S <: AbstractEquation}
    return [get_imex(equation) for equation in system]
end
##
tmp_expr = Explicit(-∂x(u*u), AB1()) + Implicit(ν * ∂x(∂x(u)) , BDF1()) + Explicit(-∂x(u*u), AB1()) + Explicit(-∂x(u*u), AB1()) 
tmp_implicit = get_implicit(tmp_expr)
tmp_implicit = get_explicit(tmp_expr)
new_eq  = get_imex(Implicit(∂t(u), BDF1()) == Explicit(-∂x(u*u), AB1()) + Implicit(ν * ∂x(∂x(u)) , BDF1()))
##
pde_system = [
    Implicit(σ, BDF1()) == Implicit(∂x(u), BDF1()),
    Implicit(∂t(u), BDF1()) == Explicit(- ∂x( u * u), AB1())+ Implicit(ν *  ∂x(σ), BDF1())
]
pde_system[1]
pde_system[2]
new_system = get_imex(pde_system)
eq1 = pde_system[1]
eq2 = pde_system[2]
println("The system of equations are")
for equation in pde_system
    println(equation)
end
println("The new system with all implicit terms 
         on the left and explicit on the right is")
for equation in new_system
    println(equation)
end