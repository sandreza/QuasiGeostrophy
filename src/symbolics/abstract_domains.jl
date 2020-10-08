
import LinearAlgebra: ×

export AbstractDomain
export AbstractBoundary

export DomainBoundary
export IntervalDomain, ProductDomain, Circle, S¹, Interval
export ×, dim, info, check_full_periodicity

abstract type AbstractDomain end
abstract type AbstractBoundary end

struct DomainBoundary <: AbstractBoundary
    closure
end
struct IntervalDomain{AT, BT, PT} <: AbstractDomain
    a::AT
    b::BT
    periodic::PT
end

function IntervalDomain(a, b; periodic=false)
    @assert a < b
    return IntervalDomain(a, b, periodic)
end
function Circle(a, b)
    @assert a < b
    return IntervalDomain(a, b, periodic = true)
end
S¹ = Circle
function Interval(a, b)
    @assert a < b
    return IntervalDomain(a, b)
end

function Base.show(io::IO, Ω::IntervalDomain) 
    a = Ω.a
    b = Ω.b
    printstyled(io, "[", color = 226)
    printstyled("$a, $b", color = 7)
    Ω.periodic ? printstyled(io, ")", color = 226) : printstyled(io, "]", color = 226)
 end


∂(Ω::IntervalDomain) = Ω.periodic ? DomainBoundary(nothing) : DomainBoundary((Ω.a, Ω.b))

struct ProductDomain{DT} <: AbstractDomain
    domains::DT
end

function Base.show(io::IO, Ω::ProductDomain) 
    for (i,domain) in enumerate(Ω.domains)
        print(domain)
        if i != length(Ω.domains)
            printstyled(io, "×", color = 118)
        end
    end
 end

 function dim(Ω::IntervalDomain)
    return 1
end

function dim(Ω::ProductDomain)
    return length(Ω.domains)
end

×(arg1::AbstractDomain, arg2::AbstractDomain) = ProductDomain((arg1, arg2))
×(args::AbstractDomain) = ProductDomain(args...)

function info(Ω::ProductDomain)
    println("This is a ", dim(Ω),"-dimensional tensor product domain.")
    print("The domain is ")
    println(Ω, ".")
    for (i,domain) in enumerate(Ω.domains)
        domain_string = domain.periodic ? "periodic" : "wall-bounded"
        println("The dimension $i domain is ", domain_string, " with length ", domain.b-domain.a)
    end
    return nothing
end

function check_full_periodicity(Ω::ProductDomain)
    b = [Ω.domains[i].periodic for i in eachindex(Ω.domains)]
    return prod(b)
end