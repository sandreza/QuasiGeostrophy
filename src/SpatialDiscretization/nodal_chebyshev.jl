"""
# Description
Julia version of Spectral Methods in Matlab

# Argument
- 'N': polynomial order

# Return
- 'D': Chebyshev differentiation matrix
- 'x': Guass-Lobatto points
"""
function cheb(N)
    if N==0
        return [0], [0]
    else
        x = @. cos(pi*(0:N)/N);
        c = [2; ones(N-1); 2] .* (-1).^(0:N);
        dX = x .- x';
        #off-diagonal entries
        D =  (c ./ c') ./ (dX + I);
        # diagonal entries
        D = D - Diagonal(sum(D', dims = 1)[1:N+1])
        return D, x
    end
end
