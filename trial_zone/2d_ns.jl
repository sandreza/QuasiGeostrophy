using QuasiGeostrophy, BenchmarkTools

Ω = S¹(0, 2π) × S¹(0, 2π)
Nx = Ny = 2^10; 
fourier_grid = FourierGrid((Nx, Ny), Ω)
fieldnames = ("u", "v", "ψ")
u, v, ψ = create_fields(names = fieldnames, grid = fourier_grid)
∂x, ∂y = create_operators(fourier_grid)
# initialize FourierField data
x, y = fourier_grid.grid
ψ(sin.(x) .* sin.(y))
u.data .= ∂y(ψ).data
v.data .= -∂x(ψ).data
Δ = FourierOperator((∂x * ∂x + ∂y * ∂y).op, FourierOperatorMetaData("Δ"))
Δ⁻¹ = inv(Δ)
Re = 0.1

function nsrhs(u,v, ∂x, ∂y, Δ, Δ⁻¹, Re)
    rhsu = -(∂x(u*u) + ∂y(u*v)) + (1/Re)*Δ(u)
    rhsv = -(∂x(u*v) + ∂y(v*v)) + (1/Re)*Δ(v)
    p = Δ⁻¹(∂x(rhsu) + ∂y(rhsv))
    rhsu -= ∂x(p)
    rhsv -= ∂y(p)
    return rhsu, rhsv
end

function nsrhs2(u,v, ∂x, ∂y, Δ, Δ⁻¹, Re)
    uu = u*u
    uv = u*v
    vu = v*u
    vv = v*v
    uux = ∂x(uu)
    uvy = ∂y(uv)
    vux = ∂x(vu)
    vvy = ∂y(vv)  
    Δu = (1/Re)*Δ(u)
    Δv = (1/Re)*Δ(v)
    rhsu = -(uux + uvy) + Δu
    rhsv = -(vux + vvy) + Δv
    p = Δ⁻¹(∂x(rhsu) + ∂y(rhsv))
    rhsu -= ∂x(p)
    rhsv -= ∂y(p)
    return rhsu, rhsv
end

P = plan_fft(u.data, flags=FFTW.PATIENT)
iP = plan_ifft(u.data, flags=FFTW.PATIENT)
function nsrhs3(u,v, ∂x, ∂y, Δ, Δ⁻¹, Re, P, iP)
    ur = iP * u.data
    vr = iP * v.data
    fmd = FourierMetaData(" ", u.metadata.grid, u.metadata.transform)
    uu = FourierField( P * (ur .* ur), fmd)
    uv = FourierField( P * (ur .* vr),fmd)
    vv = FourierField( P * (vr .* vr),fmd)
    uux = ∂x(uu)
    uvy = ∂y(uv)
    vux = ∂x(uv)
    vvy = ∂y(vv)  
    Δu = (1/Re)*Δ(u)
    Δv = (1/Re)*Δ(v)
    rhsu = -(uux + uvy) + Δu
    rhsv = -(vux + vvy) + Δv
    p = Δ⁻¹(∂x(rhsu) + ∂y(rhsv))
    rhsu -= ∂x(p)
    rhsv -= ∂y(p)
    return rhsu, rhsv
end

#nonfft operations
function srhs(u,v, ∂x, ∂y, Δ, Δ⁻¹, Re)
    uu = u
    Δu = (1/Re)*Δ(u)
    Δv = (1/Re)*Δ(v)
    rhsu = -(∂x(uu) + ∂y(uu)) + Δu
    rhsv = -(∂x(uu) + ∂y(uu)) + Δv
    p = Δ⁻¹(∂x(rhsu) + ∂y(rhsv))
    rhsu -= ∂x(p)
    rhsv -= ∂y(p)
    return rhsu, rhsv
end
