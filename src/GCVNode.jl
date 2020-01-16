module GCVNode

using ForneyLab
import ForneyLab: outboundType, isApplicable, generateId, unsafeMean, unsafeVar, @ensureVariables, addNode!, associate!, slug, format, FactorNode, SoftFactor, Interface, VariateType, unsafeMeanCov

include("./helpers.jl")

include("HGF/HGF.jl")
include("nodes/exponential_linear_quadratic.jl")
include("nodes/gaussian_controlled_variance.jl")

include("rules/definition/equality.jl")
include("rules/definition/gaussian_controlled_variance.jl")
include("rules/definition/gaussian_mean_precision.jl")

include("rules/implementation/equality.jl")
include("rules/implementation/gaussian_controlled_variance.jl")
include("rules/implementation/gaussian_mean_precision.jl")

end # module
