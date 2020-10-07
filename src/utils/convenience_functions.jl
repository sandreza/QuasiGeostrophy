export create_fields, create_operators

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
function create_fields(mod; names = (), 
                         grid = nothing, 
                         arraytype = Array,
                         floattype = ComplexF64)
    transform = Transform(grid)
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
        Base.eval(mod, :($parsed_name = FourierField($field_data, $fmd)) )
    end
    return nothing
end

"""
function create_operators(g::FourierGrid;
                          names = ("x","y","z"),
                          arraytype = Array)
# Description
A convenience function to automatically put operators in global scope

# Arguments
- `g`: FourierGrid. A fourier grid object

# Keyword Arguments
- `names`: name for operators assumes at most 3 dimensional

# Return
Nothing. But puts several Fourier Operators into the global scope
"""
function create_operators(mod, g::FourierGrid;
                          names = ("x","y","z"),
                          arraytype = Array)
    printstyled("Warning !!!", color = :red)
    println("")
    for i in 1:length(g.grid)
        operator_name = Char(0x02202) * names[i]
        print("Introducing ")
        printstyled(operator_name, color = :blue)
        println(" into the global scope.")
        k = g.wavenumbers[i]
        parsed_name = Meta.parse(operator_name)
        op = arraytype(im .* k)
        fmd = FourierOperatorMetaData(operator_name)
        Base.eval(mod, :($parsed_name = FourierOperator($op, $fmd)))
    end
    println(" ")
    return nothing
end
