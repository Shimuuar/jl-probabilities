module Roche_inplace

using Meshes
using Roots
using Interpolations
using ForwardDiff: Dual

export 
    Ω_potential,
    LagrangePoint_X, 
    Ω_critical,
    roche_r,
    StretchToRocheLobe,
    make_roche_mesh,
    InterpolatedRoche,
    apply_radial_function!,
    average_over_faces!,
    integrate_data_over_mesh,
    integrate_data_over_triangular_mesh,
    integrate_data_over_triangular_mesh2


function Ω_potential(r; mass_quotient, point_on_unit_sphere::Point)
    λ, μ, ν = coordinates(point_on_unit_sphere)
    return 1/r + 
        mass_quotient * (1 / √(1 + r^2 - 2r*λ) - λ*r) +
        (1. + mass_quotient) / 2 * r^2 * (1. - ν^2)
end

function LagrangePoint_X(mass_quotient)
    function Ω_x_derivative(x)
        -1 / x^2 - 
        mass_quotient * (1. - 1 / (x-1)^2) + 
        (1 + mass_quotient) * x
    end
    return find_zero(Ω_x_derivative, (0., 1.))
end

function Ω_critical(lagrange1_x, mass_quotient)
    return Ω_potential(lagrange1_x; mass_quotient, point_on_unit_sphere=Point(1., 0., 0.))
end

function roche_r(Ω0, lagrange1_x, mass_quotient, point_on_unit_sphere::Point)
    function Ω_partially_applied(r)
        return Ω_potential(r; mass_quotient, point_on_unit_sphere) - Ω0
    end
    find_zero(Ω_partially_applied, (0., lagrange1_x))
end




# Functions related to interpolation

struct InterpolatedRoche
    spherical_mesh
    mass_quotient_knots
    interpolants
    mesh_data
    dual_mesh_data
end


function InterpolatedRoche(spherical_mesh::SimpleMesh, mass_quotient_knots, field_names)

    r_values = zeros(nvertices(spherical_mesh), length(mass_quotient_knots))

    for (mass_quotient, r_values_for_quotient) ∈ 
            zip(mass_quotient_knots, eachcol(r_values))

        lagrange1_x = LagrangePoint_X(mass_quotient)
        Ω0 = Ω_critical(lagrange1_x, mass_quotient)
        r_values_for_quotient .= roche_r.(Ω0, lagrange1_x, mass_quotient, vertices(spherical_mesh))
    end

    interpolants = [
        cubic_spline_interpolation(mass_quotient_knots, r_values_for_vertex)
        for r_values_for_vertex ∈ eachrow(r_values)
    ]

    arrays_for_face_values = NamedTuple(field_name => zeros(nfaces(spherical_mesh, 2)) for field_name ∈ field_names)
    arrays_for_vertex_values = NamedTuple(field_name => zeros(nvertices(spherical_mesh)) for field_name ∈ field_names)
    arrays_for_vertex_values = merge(arrays_for_vertex_values, (; r = zeros(nvertices(spherical_mesh))))

    mesh_data = meshdata(
        copy(vertices(spherical_mesh)),
        topology(spherical_mesh),
        Dict(0 => arrays_for_vertex_values,
             2 => arrays_for_face_values,
        )
    )


    dual_arrays_for_face_values = NamedTuple(field_name => zeros(Dual, nfaces(spherical_mesh, 2)) for field_name ∈ field_names)
    dual_arrays_for_vertex_values = NamedTuple(field_name => zeros(Dual, nvertices(spherical_mesh)) for field_name ∈ field_names)
    dual_arrays_for_vertex_values = merge(dual_arrays_for_vertex_values, (; r = zeros(Dual, nvertices(spherical_mesh))))

    dual_mesh_data = meshdata(
        # Vector{Point{3, Dual}}(vertices(spherical_mesh)),
        Vector{Point{3, Dual}}(undef, nvertices(spherical_mesh)),
        topology(spherical_mesh),
        Dict(0 => dual_arrays_for_vertex_values,
             2 => dual_arrays_for_face_values,
        )
    )

    return InterpolatedRoche(spherical_mesh, mass_quotient_knots, interpolants, mesh_data, dual_mesh_data)
end



function stretch_mesh(mesh_data::MeshData, mass_quotient, interpolated_roche)
    r_list = values(mesh_data, 0).r
    r_list .= [
        interpolant(mass_quotient)
        for interpolant ∈ interpolated_roche.interpolants
    ]
    spherical_points = vertices(interpolated_roche.spherical_mesh)
    points = vertices(domain(mesh_data))
    @. points = Point(coordinates(spherical_points) .* r_list)
    return mesh_data
end

(interpolated_roche::InterpolatedRoche)(mass_quotient::Float64) = stretch_mesh(interpolated_roche.mesh_data, mass_quotient, interpolated_roche)

(interpolated_roche::InterpolatedRoche)(mass_quotient::Dual) = stretch_mesh(interpolated_roche.dual_mesh_data, mass_quotient, interpolated_roche)



# Functions related to filling mesh with data

function apply_radial_function!(mesh_data::MeshData, f, field_name)
    r_list = getfield(values(mesh_data, 0), :r)
    f_list = getfield(values(mesh_data, 0), field_name)
    @. f_list = f(r_list)
    return mesh_data
end

function average_over_faces!(mesh_data::MeshData, field_name)
    values_for_vertices = getfield(values(mesh_data, 0), field_name)
    values_for_faces = getfield(values(mesh_data, 2), field_name)

    values_for_faces .= map(faces(topology(domain(mesh_data)), 2)) do face
        index = indices(face)
        sum(values_for_vertices[i] for i ∈ index) / length(index)
    end

    # for (i, face) in enumerate(faces(topology(domain(mesh_data)), 2))
    #     values_for_faces[i] = sum(values_for_vertices[j] for j ∈ index) / length(index) # one(eltype(values_for_faces))
    # end

    # values_for_faces .= (
    #     avg_over_face(mesh_data, field_name, connec) 
    #     for connec in faces(topology(domain(mesh_data)), 2)
    # )

    return mesh_data
end


# Functions related to integration

function integrate_data_over_triangular_mesh(mesh_data::MeshData, field_name, direction)

    val_list = getfield(values(mesh_data, 2), field_name)
    faces_list = faces(domain(mesh_data), 2)

    return sum(val_times_area(face, val, direction) 
        for (val, face) ∈ zip(val_list, faces_list)
    )
end

function val_times_area(face, val, direction)
    # @assert isa(face, Triangle) "Для нетреугольных сеток нужно использовать integrate_data_over_mesh"
    a, b, c = vertices(face)
    n = (b-a) × (c-a)
    # s = sign(n ⋅ coordinates(a))
    dotp = n ⋅ direction # * s
    if dotp > 0
        return val * dotp / 2
    else
        return zero(val)
    end
end


function integrate_data_over_triangular_mesh2(mesh_data::MeshData, field_name, direction)
    vertices_values = getfield(values(mesh_data, 0), field_name)

    sum(faces(topology(domain(mesh_data)), 2)) do connection
        face = materialize(connection, vertices(domain(mesh_data)))
        visible_area = signed_visible_area(face, direction)
        if visible_area < 0
            return zero(visible_area)
        else
            val = avg_over_face(vertices_values, connection)
            return visible_area * val
        end
    end
end

function signed_visible_area(face, direction)
    a, b, c = vertices(face)
    n = (b-a) × (c-a)
    return (n ⋅ direction) / 2
end

function avg_over_face(vertices_values, connection)
    index = indices(connection)
    return sum(vertices_values[i] for i ∈ index) / length(index)
end


# function integrate_data_over_mesh(mesh_data::MeshData, field_name, direction)
#     val_list = getfield(values(mesh_data, 2), field_name)

#     total = 0.
#     for (val, face) ∈ zip(val_list, faces(domain(mesh_data), 2))
#         n = normal(face)
#         radius_vector = coordinates(first(vertices(face)))
#         dotp = n ⋅ direction * sign(n ⋅ radius_vector)
#         if dotp > 0
#             total += dotp * val * area(face)
#         end
#     end
#     return total
# end










# Для создания сетки один раз

# Similar to https://github.com/JuliaGeometry/Meshes.jl/blob/a0487c6824d6ee9d7389edc25ae937f1e4cf26fd/src/transforms/translate.jl

struct StretchToRocheLobe <: StatelessGeometricTransform
    mass_quotient
    lagrange1_x
    Ω0
end

function StretchToRocheLobe(mass_quotient)
    lagrange1_x = LagrangePoint_X(mass_quotient)
    Ω0 = Ω_critical(lagrange1_x, mass_quotient)
    return StretchToRocheLobe(mass_quotient, lagrange1_x, Ω0)
end

Meshes.preprocess(transform::StretchToRocheLobe, object) = transform

function Meshes.applypoint(::StretchToRocheLobe, points, prep::StretchToRocheLobe)

    function transform_point(point::Point)
        r = roche_r(prep.Ω0, prep.lagrange1_x, prep.mass_quotient, point)
        Point(coordinates(point) .* r)
    end

    return map(transform_point, points), prep
end

function make_roche_mesh(mass_quotient, discretization_method = RegularDiscretization(10))
    sphere = Sphere((0. ,0., 0.), 1.)
    return discretize(sphere, discretization_method) |>
            Rotate(Vec(0., 0., 1.), Vec(1., 0., 0.)) |>
            simplexify |>
            StretchToRocheLobe(mass_quotient)
end

end