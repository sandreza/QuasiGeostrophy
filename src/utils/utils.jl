"""
appropriate_dims(n1, n2, N)
# Description
Create an array of size n1, with the value N at the location n2, and ones elswhere

# Arguments
- `n1`: int | size of the array
- `n2`: int | location of modification
- `N` : int | value at the modification index n2

# Return
A tuple with all 1's except at location n2, where it is N
"""
function appropriate_dims(n1, n2, N)
    y = ones(Int, n1)
    y[n2] = N
    return Tuple(y)
end