using QuasiGeostrophy, LinearAlgebra, Plots, FFTW, BenchmarkTools

N = 2^12
a = 0
b = 2π
x = fourier_nodes(N, a = a, b = b)
k = fourier_wavenumbers(N, L = b - a)
f = @. sin(x) .+ 0im
plot(x, real.(f))
f̂ = copy(f)
##
FFTW.set_num_threads(1)
P = plan_fft(f)
iP = plan_ifft(f̂)
@btime mul!(f̂, P, f)
@btime mul!(f, P, f̂)
##
Nx = 2^8
Ny = 2^8
ax = 0; bx = 2π; ay = 0; by = 4π
gx  = fourier_nodes(Nx, a = ax, b = bx)
gy  = fourier_nodes(Ny, a = ay, b = by)
kx  = fourier_wavenumbers(Nx, L = bx - ax)
ky  = fourier_wavenumbers(Ny, L = by - ay)
x = reshape(gx, (Nx, 1))
y = reshape(gy, (1, Ny))
kx = reshape(kx, (Nx, 1))
ky = reshape(ky, (1, Ny))
f = @. 2 * sin(x) + sin(2*y) + 0im
∂ˣf = @. 2 * cos(x) + 0im * y
∂ʸf = @. 2 * cos(2*y) + 0im * x
Δf = @. - 2 * sin(x) - 4 * sin(2*y)
f̂ = copy(f)
FFTW.set_num_threads(Threads.nthreads())
P = plan_fft(f)
# typeof(iP) <: AbstractFFTs.ScaledPlan
iP = plan_ifft(f)
mul!(f̂, P, f)
# @btime mul!(f̂, P, f)
##
contourf(gy, gx, real.(f) )

##
import Base: *, ^, +, inv
abstract type Filter end
abstract type NoFilter <: Filter end
abstract type OrszagFilter <: Filter end
struct PartialDerivative{T}
    op::T
end
function *(∂x::PartialDerivative, ϕ::AbstractArray)
    return ∂x.op .* ϕ
end
function (p::PartialDerivative)(ϕ::AbstractArray)
    return *(p, ϕ)
end
function ^(p::PartialDerivative, α::Number)
    return PartialDerivative(p.op.^(α))
end
function +(p::PartialDerivative, q::PartialDerivative)
    return PartialDerivative(p.op .+ q.op)
end
function convolve(û, v̂, P, iP)
    u = iP * û
    v = iP * v̂
    return P * (u .* v)
end
function box_filter(û)
    n = length(û)
    mid = div(n, 2) + 1
    pm  = div(n, 6)
    u = copy(û)
    u[mid - pm : mid + pm] .= 0.0
    return u
end
function filter_convolve(û, v̂, P, iP)
    u = box_filter(û)
    v = box_filter(v̂)
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
function inv(a::PartialDerivative)
    inv_op = 1 ./ a.op 
    @inbounds for i in eachindex(inv_op)
        if abs(inv_op[i]) == Inf
            inv_op[i] = 1.0
        elseif isnan(norm(inv_op[i]))
            inv_op[i] = 1.0
        end
    end
    return PartialDerivative(inv_op)
end
##
∂x = PartialDerivative(im .* kx)
∂y = PartialDerivative(im .* ky)
Δ = ∂x^2 + ∂y^2

norm(∂ˣf - iP * ∂x(f̂)) ./ norm(∂ˣf )
norm(∂ʸf - iP * ∂y(f̂)) ./ norm(∂ʸf )
norm(Δf  - iP * Δ(f̂) ) ./ norm(Δf)
