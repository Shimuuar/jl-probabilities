module LuminocityModels

using Turing
using StructTypes
using Revise

include("Roche.jl")
using .Roche

export
    MeshParams,
    ModelParams,
    ChainParams,
    first_model,
    star_magnitude,
    T_4,
    planck_formula,
    black_body_K_rectangle


@kwdef struct MeshParams
    n_discretization_points::Int = 64
    mass_quotient_nodes = 0.1 : 0.1 : 10.
end

@kwdef struct ModelParams
    mesh_params::MeshParams = MeshParams()
    model_function::Function = first_model
    period::Union{Float64, Nothing}
    β::Union{Nothing, Float64}
    fixed_σ::Union{Nothing, Float64, Vector{Float64}} = nothing
    luminocity_function::Function = T_4
    fixed_temperature_at_bottom::Union{Nothing, Float64} = nothing
    measurements_t::Vector{Float64}
    measurements_y::Vector{Float64}
end

@kwdef struct ChainParams
    model_params::ModelParams
    n_samples::Int
    n_chains::Int = 1
    sampler::Turing.InferenceAlgorithm = NUTS()
    init_params::Union{Nothing, NamedTuple} = nothing
end
StructTypes.StructType(::Turing.InferenceAlgorithm) = StructTypes.StringType()


T_4(T) = T^4

function planck_formula(λ, T)
    h = 6.62607004e-34
    c = 299792458
    k = 1.38064852e-23
    return 2h*c^2 / (λ^5 * (exp(h*c / (λ*k*T)) - 1))
end

function black_body_K_rectangle(T)
    λ = 2.2e-6
    Δλ = 0.4e-6
    return planck_formula(λ, T) * Δλ
end


StructTypes.StructType(::typeof(T_4)) = StructTypes.StringType()
StructTypes.StructType(::typeof(black_body_K_rectangle)) = StructTypes.StringType()


function first_model(model_params)
    interpolated_mesh = InterpolatedRocheMesh(
        model_params.mesh_params.n_discretization_points,
        model_params.mesh_params.mass_quotient_nodes
    )

    mass_quotient_min = model_params.mesh_params.mass_quotient_nodes[1]
    mass_quotient_max = model_params.mesh_params.mass_quotient_nodes[end]


    @model function model(measurements_t, measurements_y)
        initial_phase ~ Uniform(-π, π)
        phases = initial_phase .+ measurements_t * 2π / model_params.period

        mass_quotient ~ Uniform(mass_quotient_min, mass_quotient_max)
        observer_angle ~ Uniform(0., π)

        if model_params.fixed_temperature_at_bottom === nothing
            temperature_at_bottom ~ FlatPos(500.)
        else
            temperature_at_bottom = model_params.fixed_temperature_at_bottom
        end

        predicted_magnitudes = star_magnitude(
            phases;
            mass_quotient,
            observer_angle,
            temperature_at_bottom,
            model_params.β,
            interpolated_mesh,
            model_params.luminocity_function
        )

        offset ~ Flat()
        predicted_magnitudes .+= offset

        if model_params.fixed_σ === nothing
            σ ~ FlatPos(0.)
        else
            σ = model_params.fixed_σ
        end

        measurements_y .~ Normal.(predicted_magnitudes, σ)
    end

    return model(model_params.measurements_t, model_params.measurements_y)
end

StructTypes.StructType(::typeof(first_model)) = StructTypes.StringType()



function star_magnitude(phases; mass_quotient, observer_angle,
                        temperature_at_bottom, β, interpolated_mesh,
                        luminocity_function)

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

    luminocities = [
        integrate_data_over_triangular_mesh(mesh, :L, direction)
        for direction ∈ directions
    ]
    return @. -2.5 * log10(luminocities)
end

end