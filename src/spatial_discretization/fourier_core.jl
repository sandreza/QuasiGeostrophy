include("fourier_utils.jl")
include("fourier_operations.jl")
include("fourier_fields.jl")

# Define Interactions
function (∇::FourierOperator)(ϕ::FourierField)
    return FourierField(∇.op .* ϕ.data, ϕ.metadata)
end
# Define Naming
function (∇::FourierOperator{S, T})(ϕ::FourierField) where
    {S, T <: FourierOperatorMetaData}
    op_name = ∇.metadata.name
    ϕ_name = ϕ.metadata.name
    new_name = op_name * "(" * ϕ_name * ")"
    fmd = FourierMetaData(new_name, ϕ.metadata.grid, ϕ.metadata.transform)
    return FourierField(∇.op .* ϕ.data, fmd)
end