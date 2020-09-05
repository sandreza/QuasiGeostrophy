using FFTW

"""
fourier_points(n; a = 0, b = 2π)

# Description
- Create a uniform grid of points for periodic functions

# Arguments
- `N`: integer | number of evenly spaced points 

# Keyword Arguments
- `a`: number | starting point of interval [a, b)
- `b`: number | ending point of interval [a, b)

# Return
- `g`: array | an array of points of evenly spaced points from [a, b)
"""
function fourier_nodes(N; a = 0, b = 2π)
    return (b-a) .* collect(0:(N-1))/N .+ a
end


"""
wavevec(N; L = 2π)

# Description
- Create wavenumbers associated with the domain of length L

# Arguments
- `N`: integer | number of wavevectors

# Keyword Arguments
- `L`: number | length of interval [a, b), L = b-a

# Return
- `wavenumbers`: array | an array of wavevectors
"""
function fourier_wavenumbers(N; L = 2π)
    up = collect(0:1:N-1)
    down = collect(-N:1:-1)
    indices = up
    indices[div(N,2)+1:end] = down[div(N,2)+1:end]
    wavenumbers = 2π / L .* indices
    return wavenumbers
end