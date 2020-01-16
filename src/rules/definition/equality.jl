mutable struct SPEqualityGaussianGCV <: SumProductRule{Equality} end
outboundType(::Type{SPEqualityGaussianGCV}) = Message{GaussianMeanVariance}
function isApplicable(::Type{SPEqualityGaussianGCV}, input_types::Vector{Type})
    nothing_inputs = 0
    gaussian_inputs = 0
    exp_lin_quad_inputs = 0
    for input_type in input_types
        if input_type == Nothing
            nothing_inputs += 1
        elseif matches(input_type, Message{Gaussian})
            gaussian_inputs += 1
        elseif matches(input_type, Message{ExponentialLinearQuadratic})
            exp_lin_quad_inputs += 1
        end
    end

    return (nothing_inputs == 1) && (gaussian_inputs == 1) && (exp_lin_quad_inputs == 1)
end
