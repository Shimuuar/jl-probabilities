module LuminocityModels

using Turing
using StructTypes
using Meshes
using Revise

include("Roche.jl")
include("LuminocityFunctions.jl")

using .Roche
using .LuminocityFunctions

export
    MeshParams,
    ChannelParams,
    ModelParams,
    ChainParams,
    first_model,
    star_magnitude


@kwdef struct MeshParams
    n_discretization_points::Int = 64
    mass_quotient_nodes = 0.1 : 0.1 : 10.
end

@kwdef struct ChannelParams
    measurements_t::Vector{Float64}
    measurements_y::Vector{Float64}
    darkening_function = one1
    darkening_coefficients = ()
    luminocity_function = T_4
    σ_measured::Vector{Float64}
    σ_common::Union{Float64, Distribution} = 0.
end

@kwdef struct ModelParams
    mesh_params::MeshParams = MeshParams()
    model_function::Function = first_model
    channels::Vector{ChannelParams}
    period::Union{Float64, Distribution}
    β::Union{Float64, Distribution} = 0.25
    temperature_at_bottom::Union{Float64, Distribution} = 3500.
end
StructTypes.StructType(::Distribution) = StructTypes.StringType()

@kwdef struct ChainParams
    model_params::ModelParams
    n_samples::Int
    n_chains::Int = 1
    sampler::Turing.InferenceAlgorithm = NUTS()
    init_params::Union{Nothing, NamedTuple} = nothing
end
StructTypes.StructType(::Turing.InferenceAlgorithm) = StructTypes.StringType()



macro dist_or_const(var, val)
    esc(:(
        if isa($val, Distribution)
            $var ~ $val
        elseif isa($val, Float64) || isa($val, Vector{Float64})
            $var = $val
        else
            throw(ArgumentError("val must be Distribution, Float64 or Vector{Float64}"))
        end
    ))
end


function first_model(model_params)
    interpolated_mesh = InterpolatedRocheMesh(
        model_params.mesh_params.n_discretization_points,
        model_params.mesh_params.mass_quotient_nodes
    )

    mass_quotient_min = model_params.mesh_params.mass_quotient_nodes[1]
    mass_quotient_max = model_params.mesh_params.mass_quotient_nodes[end]


    @model function model(channels::Vector{ChannelParams}, measurements_y = Float64[])
        mass_quotient ~ Uniform(mass_quotient_min, mass_quotient_max)
        observer_angle ~ Uniform(0., π/2)

        @dist_or_const temperature_at_bottom model_params.temperature_at_bottom
        @dist_or_const β model_params.β

        initial_phase ~ Uniform(-π, π)

        @dist_or_const period model_params.period

        offset ~ filldist(Flat(), length(channels))
        σ_common = Array{Float64}(undef, length(channels))

        for (i, channel) ∈ enumerate(channels)
            phases = initial_phase .+ channel.measurements_t * 2π / period

            predicted_magnitudes = star_magnitude(
                phases;
                mass_quotient,
                observer_angle,
                temperature_at_bottom,
                β,
                interpolated_mesh,
                channel.luminocity_function,
                channel.darkening_function,
                channel.darkening_coefficients
            )

            predicted_magnitudes .+= offset[i]

            @dist_or_const σ_common[i] channel.σ_common
            σ = @. √(channel.σ_measured^2 + σ_common[i]^2)

            measurements_y = channel.measurements_y
            measurements_y .~ Normal.(predicted_magnitudes, σ)
        end
    end

    return model(model_params.channels)
end

StructTypes.StructType(::typeof(first_model)) = StructTypes.StringType()


function star_magnitude(phases; mass_quotient, observer_angle,
                        temperature_at_bottom, β, interpolated_mesh,
                        luminocity_function, darkening_function, darkening_coefficients)

    directions = [(
        sin(observer_angle) * cos(phase),
        sin(observer_angle) * sin(phase),
        cos(observer_angle)
    ) for phase ∈ phases]

    function temperature(g)
        return temperature_at_bottom * abs(g)^β
    end

    mesh = interpolated_mesh(mass_quotient)
    mesh = apply_function(mesh, temperature, :g, :T)
    mesh = apply_function(mesh, luminocity_function, :T, :L)

    normals = calc_function_on_faces(mesh, normal)
    areas = calc_function_on_faces(mesh, area)

    luminocities = [
        integrate_data_over_mesh(mesh, :L, direction, normals, areas, darkening_function, darkening_coefficients)
        for direction ∈ directions
    ]
    return @. -2.5 * log10(luminocities)
end

end