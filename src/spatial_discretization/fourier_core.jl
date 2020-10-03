include("fourier_utils.jl")
include("fourier_operations.jl")
include("fourier_fields.jl")

# Define Interactions
function (∇::FourierOperator)(ϕ::FourierField)
    return FourierField(∇.op .* ϕ.data, ϕ.metadata)
end