
function Operator(O::FourierOperator{S, T}) where
    {S, T <: FourierOperatorMetaData}
    omd = OperatorMetaData(O, O.metadata.name)
    return Operator(nothing, omd)
end

function compute(O::Operator{S, OperatorMetaData{T, U}}) where
    {S <: Nothing, T <: FourierOperator, U}
    return O.metadata.operation
end

