# Quick Structs for checking calculations
using QuasiGeostrophy
import QuasiGeostrophy:compute

export Wrapper, @wrapper
export DirectionalDerivative, GradientMetaData

struct Wrapper{T, S} <: AbstractExpression
    data::T
    meta_data::S
end
# Struct for MetaData
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
##

struct DirectionalDerivative{𝒮} <: AbstractExpression
    direction::𝒮
end
struct GradientMetaData{𝒮}
    direction::𝒮
end
function (p::DirectionalDerivative)(expr::AbstractExpression)
    return Gradient(expr, GradientMetaData(p.direction))
end
function Base.show(io::IO, p::DirectionalDerivative{S}) where S <: String
    print(io, Char(0x02202) * p.direction)
end
function Base.show(io::IO, p::Gradient{S, T}) where {S, T <: GradientMetaData{String}}
    printstyled(io, Char(0x02202) * p.metadata.direction, "(", color = 165)
    print(io, p.operand)
    printstyled(io, ")", color = 165)
end

∂x = DirectionalDerivative("x")
∂y = DirectionalDerivative("y")
∂z = DirectionalDerivative("z")
∂t = DirectionalDerivative("t")

## Probably should move to main directory
size(f::FourierGrid) = length.(f.grid)

"""
function create_fields(; names = (), 
                         grid = nothing, 
                         transform = nothing,
                         arraytype = Array,
                         floattype = ComplexF64)
# Description
Automates the constructions of fourier fields with names

# Arguments
All the arguments are keyword arguments

# Keyword Arguments
- `names`: tuple of strings
- `grid`: FourierGrid object
- `transfrom`: Transform object
- `arraytype`: default = Array, 
- `floattype`: default = ComplexF64

# Return
nothing
"""
function create_fields(; names = (), 
                         grid = nothing, 
                         transform = nothing,
                         arraytype = Array,
                         floattype = ComplexF64)
    printstyled("Warning!!! ", color = :red)
    print("the name(s) ")
    printstyled("(", color = :blue)
    for name in names 
        printstyled(name, ", ", color = :blue)
    end
    printstyled(")", color = :blue)
    print(" are being overwritten in the ")
    print("global scope.")
    println(" ")
    dimensions = size(grid)
    for name in names
        local fmd = FourierMetaData(name, grid, transform)
        local field_data = arraytype(zeros(floattype, dimensions...))
        local parsed_name = Meta.parse(string(name))
        @eval $parsed_name = FourierField($field_data, $fmd)
    end
    return nothing
end

# convenience function to initialize array
function (ϕ::FourierField)(a::AbstractArray)
    f1 = ϕ.metadata.transform.forward * (a .+ 0im)
    ϕ.data .= f1
    return nothing
end

import LinearAlgebra: norm
norm(ϕ::FourierField) = norm(ϕ.data)

function spectrum(ϕ::FourierField)
    f = ϕ.data
    g = log10.(abs.(f))[1:div(length(f),2)+1]
    ymax = maximum(g) * 1.1
    ymin = 1e-16 * ymax
    wvi = collect(1:div(length(f),2)+1) .- 1
    p1 = scatter(wvi, g, label = ϕ.metadata.name)
    scatter!(xlabel = "wavenumber index")
    scatter!(ylabel = "log10(spectral amplitude)")
    scatter!(ylims = (ymin, ymax))
    return p1
end

function create_operators(g::FourierGrid;
                          names = ("x","y","z"),
                          arraytype = Array)
    printstyled("Warning !!!", color = :red)
    for i in 1:length(g.grid)
        operator_name = Char(0x02202) * names[i]
        print(" introducing ")
        printstyled(operator_name, color = :blue)
        println(" into the global scope")
        k = fourier_grid.wavenumbers[i]
        parsed_name = Meta.parse(operator_name)
        op = arraytype(im .* k)
        fmd = FourierOperatorMetaData(operator_name)
        @eval $parsed_name = FourierOperator($op, $fmd)
    end
    println(" ")
    return nothing
end
