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
    ModelParams,
    ChainParams,
    first_model,
    star_magnitude


@kwdef struct MeshParams
    n_discretization_points::Int = 64
    mass_quotient_nodes = 0.1 : 0.1 : 10.
end

@kwdef struct ModelParams
    mesh_params::MeshParams = MeshParams()
    model_function::Function = first_model
    period::Union{Float64, Distribution}
    β::Union{Float64, Distribution} = 0.25
    σ::Union{Float64, Vector{Float64}, Distribution} = 0.1
    luminocity_function = T_4
    temperature_at_bottom::Union{Float64, Distribution} = 3500.
    darkening_function = one1
    darkening_coefficients = ()
    measurements_t::Vector{Float64}
    measurements_y::Vector{Float64}
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
        elseif isa($val, Float64)
            $var = $val
        else
            throw(ArgumentError("val must be either Distribution or Float64"))
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


    @model function model(measurements_t, measurements_y)
        initial_phase ~ Uniform(-π, π)

        @dist_or_const period model_params.period

        phases = initial_phase .+ measurements_t * 2π / period

        mass_quotient ~ Uniform(mass_quotient_min, mass_quotient_max)
        observer_angle ~ Uniform(0., π)

        @dist_or_const temperature_at_bottom model_params.temperature_at_bottom
        @dist_or_const β model_params.β

        predicted_magnitudes = star_magnitude(
            phases;
            mass_quotient,
            observer_angle,
            temperature_at_bottom,
            β,
            interpolated_mesh,
            model_params.luminocity_function,
            model_params.darkening_function,
            model_params.darkening_coefficients
        )

        offset ~ Flat()
        predicted_magnitudes .+= offset

        @dist_or_const σ model_params.σ

        measurements_y .~ Normal.(predicted_magnitudes, σ)
    end

    return model(model_params.measurements_t, model_params.measurements_y)
end

StructTypes.StructType(::typeof(first_model)) = StructTypes.StringType()


function second_model(model_params)
    interpolated_mesh = InterpolatedRocheMesh(
        model_params.mesh_params.n_discretization_points,
        model_params.mesh_params.mass_quotient_nodes
    )

    mass_quotient_min = model_params.mesh_params.mass_quotient_nodes[1]
    mass_quotient_max = model_params.mesh_params.mass_quotient_nodes[end]


    @model function model(measurements_t, measurements_y)
        initial_phase ~ Uniform(-π, π)

        @dist_or_const period model_params.period

        phases = initial_phase .+ measurements_t * 2π / period

        mass_quotient ~ Uniform(mass_quotient_min, mass_quotient_max)
        observer_angle ~ Uniform(0., π)

        temperature_at_bottom ~ model_params.temperature_at_bottom
        @dist_or_const β model_params.β

        predicted_magnitudes = star_magnitude(
            phases;
            mass_quotient,
            observer_angle,
            temperature_at_bottom,
            β,
            interpolated_mesh,
            model_params.luminocity_function,
            model_params.darkening_function,
            model_params.darkening_coefficients
        )

        offset ~ Flat()
        predicted_magnitudes .+= offset

        @dist_or_const σ model_params.σ

        measurements_y .~ Normal.(predicted_magnitudes, σ)
    end

    return model(model_params.measurements_t, model_params.measurements_y)
end

StructTypes.StructType(::typeof(second_model)) = StructTypes.StringType()



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