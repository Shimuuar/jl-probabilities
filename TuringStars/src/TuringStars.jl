module TuringStars

using Reexport
using Revise

include("LuminocityModels.jl")
include("Roche.jl")
include("ChainCache.jl")

@reexport using .Roche
@reexport using .LuminocityModels
@reexport using .ChainCache


function barr()
    return("Hello, world!")
end

# export
#     Ω_potential,
#     Ω_grad,
#     LagrangePoint_X,
#     Ω_critical,
#     roche_r,
#     StretchToRocheLobe,
#     make_roche_geotable,
#     InterpolatedRocheMesh,
#     integrate_data_over_triangular_mesh,
#     apply_radial_function,
#     apply_function,

#     MeshParams,
#     ModelParams,
#     ChainParams,
#     model_from_params,
#     star_magnitude,
#     T_4,
#     planck_formula,
#     black_body_K_rectangle,

#     cached_sample

end
