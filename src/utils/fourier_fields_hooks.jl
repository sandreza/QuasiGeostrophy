
for unary_operator in unary_operators
    b_symbol = Meta.parse.(unary_operator[2]) #broadcast
    @eval import Base: $b_symbol
    @eval function $b_symbol(field1::FourierField)
        data = broadcast($b_symbol, field1.data)
        metadata  = field1.metadata
        symbname = string($b_symbol)
        name = symbname * "(" * field1.metadata.name * ")"
        fmd = FourierMetaData(name, metadata.grid, metadata.transform)
        FourierField(data, fmd )
    end
end

for binary_operator in [binary_operators..., ["Negative", "-"]]
    b_symbol = Meta.parse.(binary_operator[2]) 
    @eval import Base: $b_symbol
    @eval function $b_symbol(field1::FourierField, field2::𝒮) where {𝒮}
        data = broadcast($b_symbol, field1.data, field2)
        metadata  = field1.metadata
        symbname = string($b_symbol)
        name1 = field1.metadata.name 
        name2 = string(field2)
        name = "(" * name1 * symbname * name2 * ")"
        fmd = FourierMetaData(name, metadata.grid, metadata.transform)
        return FourierField(data, fmd)
    end
    @eval function $b_symbol(field1::𝒯, field2::FourierField) where {𝒯}
        data = broadcast($b_symbol, field1, field2.data)
        metadata  = field2.metadata
        symbname = string($b_symbol)
        name1 = string(field1)
        name2 = field2.metadata.name 
        name = "(" * name1 * symbname * name2 * ")"
        fmd = FourierMetaData(name, metadata.grid, metadata.transform)
        return FourierField(data, fmd)
    end
end

# exception for multiplication of fourier fields
for binary_operator in [binary_operators[1], ["Negative", "-"]]
    b_symbol = Meta.parse.(binary_operator[2]) 
    @eval function $b_symbol(field1::FourierField, field2::FourierField)
        data = broadcast($b_symbol, field1.data, field2.data)
        metadata  = field1.metadata
        symbname = string($b_symbol)
        name1 = field1.metadata.name 
        name2 = field2.metadata.name 
        name = "(" * name1 * symbname * name2 * ")"
        fmd = FourierMetaData(name, metadata.grid, metadata.transform)
        return FourierField(data, fmd)
    end
end

function *(f̂::FourierField, ĝ::FourierField)
    fwd = f̂.metadata.transform.forward
    bwd = f̂.metadata.transform.backward
    f = real.(bwd * f̂.data)
    g = real.(bwd * ĝ.data)
    fg = broadcast(*, f, g)
    metadata  = ĝ.metadata
    name1 = f̂.metadata.name
    name2 = ĝ.metadata.name 
    name = "(" * name1 * "*" * name2 * ")"
    fmd = FourierMetaData(name, metadata.grid, metadata.transform)
    return FourierField(fwd * fg, fmd)
end
function ^(f̂::FourierField, n::Number)
    fwd = f̂.metadata.transform.forward
    bwd = f̂.metadata.transform.backward
    f = real.(bwd * f̂.data)
    ff = broadcast(^, f, n)
    metadata  = f̂.metadata
    name1 = f̂.metadata.name
    name = "(" * name1 * ")" * "^" * string(n)
    fmd = FourierMetaData(name, metadata.grid, metadata.transform)
    return FourierField(fwd * ff, fmd)
end

compute(a::FourierField) = a

## need the following hook
function compute(a::DerivativeMetaData{𝒮,𝒯}) where
    {𝒮 <: FourierOperator, 𝒯}
    return a.operation
end 