export create_fields, create_operators

"""
function create_fields(; names = (), 
                         grid = nothing, 
                         transform = nothing,
                         arraytype = Array,
                         floattype = ComplexF64)
# Description
Automates the constructions of fourier fields with names.

# Arguments
none

# Keyword Arguments
- `names`: tuple of strings
- `grid`: FourierGrid object
- `transfrom`: Transform object
- `arraytype`: default = Array, 
- `floattype`: default = ComplexF64

# Return
- Array of fourier field objects with names
"""
function create_fields(; names = (), 
                         grid = nothing, 
                         arraytype = Array,
                         floattype = ComplexF64)
    transform = Transform(grid)
    dimensions = size(grid)
    state = []
    for name in names
        fmd = FourierMetaData(name, grid, transform)
        field_data = arraytype(zeros(floattype, dimensions...))
        push!(state, FourierField(field_data, fmd))
    end
    if length(state) == 1
        return state[1]
    else
        return state
    end
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
several operators
"""
function create_operators(g::FourierGrid;
                          names = ("x","y","z"),
                          arraytype = Array)
    @assert length(names) <= 3
    ops = []
    for i in 1:length(g.grid)
        operator_name = Char(0x02202) * names[i]
        k = g.wavenumbers[i]
        op = arraytype(im .* k)
        fmd = FourierOperatorMetaData(operator_name)
        push!(ops, FourierOperator(op, fmd))
    end
    if length(ops) == 1
        return ops[1]
    else
        return ops
    end
end

## DO NOT USE THESE FUNCTIONS NORMALLY !!!!!
"""
function create_fields(mod::Module; names = (), 
                         grid = nothing, 
                         transform = nothing,
                         arraytype = Array,
                         floattype = ComplexF64)
# Description
Automates the constructions of fourier fields with names
# Arguments
- `mod`: Module for evaluating functions. Typically mod = @__MODULE__
# Keyword Arguments
- `names`: tuple of strings
- `grid`: FourierGrid object
- `transfrom`: Transform object
- `arraytype`: default = Array, 
- `floattype`: default = ComplexF64
# Return
nothing
"""
function create_fields(mod::Module; names = (), 
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
function create_operators(mod::Module, g::FourierGrid;
                          names = ("x","y","z"),
                          arraytype = Array)
# Description
A convenience function to automatically put operators in global scope.
Should only be used for quick tests.

# Arguments
- `mod`: Module for evaluating functions. Typically mod = @__MODULE__
- `g`: FourierGrid. A fourier grid object
# Keyword Arguments
- `names`: name for operators assumes at most 3 dimensional
# Return
Nothing. But puts several Fourier Operators into the global scope
"""
function create_operators(mod::Module, g::FourierGrid;
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