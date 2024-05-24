module LuminocityModels

using Turing
using StructTypes
using Meshes
using Revise

# using Metal

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
    inverted_model,
    star_magnitude,
    star_magnitude_GPU


@kwdef struct MeshParams
    catmullclark_iterations::Int = 4
    mass_quotient_nodes = 0.1 : 0.1 : 10.
end

@kwdef struct ChannelParams
    measurements_t::Vector{Float64}
    measurements_y::Vector{Float64}
    darkening_function = one1
    darkening_coefs_interpolant = T -> ()
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
        tetra_sphere(model_params.mesh_params.catmullclark_iterations),
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
                channel.darkening_coefs_interpolant
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


function inverted_model(model_params)
    interpolated_mesh = InterpolatedRocheMesh(
        tetra_sphere(model_params.mesh_params.catmullclark_iterations),
        model_params.mesh_params.mass_quotient_nodes
    )

    mass_quotient_min = 1 / model_params.mesh_params.mass_quotient_nodes[end]
    mass_quotient_max = 1 / model_params.mesh_params.mass_quotient_nodes[1]


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
                mass_quotient = 1 / mass_quotient,
                observer_angle,
                temperature_at_bottom,
                β,
                interpolated_mesh,
                channel.luminocity_function,
                channel.darkening_function,
                channel.darkening_coefs_interpolant
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

StructTypes.StructType(::typeof(inverted_model)) = StructTypes.StringType()


function star_magnitude(phases; mass_quotient, observer_angle,
                        temperature_at_bottom, β, interpolated_mesh,
                        luminocity_function, darkening_function, darkening_coefs_interpolant)

    directions = [(
        sin(observer_angle) * cos(phase),
        sin(observer_angle) * sin(phase),
        cos(observer_angle)
    ) for phase ∈ phases]

    function temperature(g)
        return temperature_at_bottom * abs(g)^β
    end

    mesh = interpolated_mesh(mass_quotient)
    mesh = avg_over_faces(mesh, :g)
    mesh = apply_function(mesh, temperature, :g, :T)
    mesh = apply_function(mesh, luminocity_function, :T, :L)
    mesh = apply_function(mesh, darkening_coefs_interpolant, :T, :darkening_coefs)

    normals = calc_function_on_faces(mesh, normalized_normal)
    areas = calc_function_on_faces(mesh, area)

    luminocities = [
        integrate_data_over_mesh(mesh, :L, direction, normals, areas, darkening_function)
        for direction ∈ directions
    ]
    return @. -2.5 * log10(luminocities)
end





# function star_magnitude_GPU(phases; mass_quotient, observer_angle,
#                         temperature_at_bottom, β, interpolated_mesh,
#                         luminocity_function, darkening_function, darkening_coefs_interpolant)

#     directions = [(
#         Float32(sin(observer_angle) * cos(phase)),
#         Float32(sin(observer_angle) * sin(phase)),
#         Float32(cos(observer_angle))
#     ) for phase ∈ phases] |> mtl

#     function temperature(g)
#         return temperature_at_bottom * abs(g)^β
#     end

#     mesh = interpolated_mesh(mass_quotient)
#     mesh = avg_over_faces(mesh, :g)
#     mesh = apply_function(mesh, temperature, :g, :T)
#     mesh = apply_function(mesh, luminocity_function, :T, :L)
#     mesh = apply_function(mesh, darkening_coefs_interpolant, :T, :darkening_coefs)

#     normals = calc_function_on_faces(mesh, normalized_normal) .|> Meshes.Vec3f |> mtl
#     areas = calc_function_on_faces(mesh, area) |> mtl
#     Ls = values(mesh, 2).L |> mtl
#     darkening_coefs = values(mesh, 2).darkening_coefs .|> Meshes.Vec{4, Float32} |> mtl

#     luminocities = Metal.zeros(length(directions), length(normals))

#     function add_face(directions, normals, areas, Ls, darkening_coefs)
#         direction_no, face_no = thread_position_in_grid_2d()
#         if direction_no > length(directions) || face_no > length(normals)
#             return
#         end
#         direction = directions[direction_no]
#         normal = normals[face_no]
#         area = areas[face_no]
#         L = Ls[face_no]
#         darkening_coef = darkening_coefs[face_no]

#         cosine = direction ⋅ normal
#         if cosine < zero(cosine)
#             return
#         else
#             darkening = darkening_function(cosine, darkening_coef...)
#             luminocities[direction_no, face_no] += cosine * area * L * darkening
#             return
#         end
#     end

#     grid = (length(directions), length(normals))
#     @metal groups=grid add_face(directions, normals, areas, Ls, darkening_coefs)

#     luminocities_ = sum(luminocities, dims=2)

#     return Float32(-2.5) .* log10.(luminocities_)
# end



end