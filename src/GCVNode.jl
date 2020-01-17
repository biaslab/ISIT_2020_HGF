module GCVNode

using ForneyLab

include("./imports.jl")

include("./helpers.jl")

# ForneyLab node definitions
include("HGF/HGF.jl")
include("nodes/exponential_linear_quadratic.jl")
include("nodes/gaussian_controlled_variance.jl")

# Update rules signature registration
include("rules/definition/equality.jl")
include("rules/definition/gaussian_controlled_variance.jl")
include("rules/definition/gaussian_mean_precision.jl")

# Update rules implementation
include("rules/implementation/equality.jl")
include("rules/implementation/gaussian_controlled_variance.jl")
include("rules/implementation/gaussian_mean_precision.jl")

end # module
