module Roche

using Meshes
using GeoTables
using Roots
using Interpolations
using LinearAlgebra

export
    Ω_potential,
    LagrangePoint_X, 
    Ω_critical,
    roche_r,
    StretchToRocheLobe,
    make_roche_geotable,
    InterpolatedRocheMesh,
    integrate_data_over_triangular_mesh,
    apply_radial_function


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
    if point_on_unit_sphere.coords[1] ≈ one(eltype(point_on_unit_sphere.coords))
        return lagrange1_x
    else
        return find_zero(Ω_partially_applied, (0., lagrange1_x))
    end
end




# Functions related to interpolation

struct InterpolatedRocheMesh
    spherical_mesh
    mass_quotient_knots
    interpolants
    r_lists
    points_lists
end


function InterpolatedRocheMesh(spherical_mesh::SimpleMesh, mass_quotient_knots)

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

    return InterpolatedRocheMesh(spherical_mesh, mass_quotient_knots, interpolants, Dict(), Dict())
end

function _get_r_list(interpolated_mesh::InterpolatedRocheMesh, ::Type{T}) :: Vector{T} where T
    if haskey(interpolated_mesh.r_lists, T)
        return interpolated_mesh.r_lists[T]
    else
        r_list = zeros(T, length(interpolated_mesh.mass_quotient_knots))
        interpolated_mesh.r_lists[T] = r_list
        return r_list
    end
end

function _get_points_list(interpolated_mesh::InterpolatedRocheMesh, ::Type{T}) :: Vector{Point{3, T}} where T
    if haskey(interpolated_mesh.points_lists, T)
        return interpolated_mesh.points_lists[T]
    else
        points_list = zeros(Point{3, T}, length(interpolated_mesh.mass_quotient_knots))
        interpolated_mesh.points_lists[T] = points_list
        return points_list
    end
end

function (interpolated_mesh::InterpolatedRocheMesh)(mass_quotient::T) where T
    points::Vector{Point3} = vertices(interpolated_mesh.spherical_mesh)

    r_list = _get_r_list(interpolated_mesh, T)
    r_list .= mass_quotient .|> interpolated_mesh.interpolants

    new_points = _get_points_list(interpolated_mesh, T)
    new_points .= [
        Point(coordinates(point) .* r)
        for (point, r) ∈ zip(points, interpolated_mesh.r_lists[T])
    ]

    new_mesh::SimpleMesh{3, T, Vector{Point{3, T}}, SimpleTopology{Connectivity{Triangle, 3}}} = SimpleMesh(
        new_points,
        topology(interpolated_mesh.spherical_mesh)
    )
    return GeoTable(
        new_mesh,
        Dict(0 => (r = r_list,))
    )
end



# Functions related to integration


function integrate_data_over_triangular_mesh(geo_table::GeoTable, field_name, direction)
    vertices_values = getfield(values(geo_table, 0), field_name)

    sum(faces(topology(domain(geo_table)), 2)) do connection
        face = materialize(connection, vertices(domain(geo_table)))
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


function apply_radial_function(geo_table::GeoTable, f, field_name)
    r_list = values(geo_table, 0).r
    f_list = f.(r_list)

    return GeoTable(
        domain(geo_table),
        Dict(0 => (r = r_list, field_name => f_list,))
    ) 
end







# Для создания сетки один раз

function make_roche_geotable(mass_quotient, discretization_method = RegularDiscretization(10))
    sphere = Sphere((0. ,0., 0.), 1.)
    spherical_mesh = discretize(sphere, discretization_method) |>
                    Rotate(Vec(0., 0., 1.), Vec(1., 0., 0.)) |>
                    simplexify
    lagrange1_x = LagrangePoint_X(mass_quotient)
    Ω0 = Ω_critical(lagrange1_x, mass_quotient)
    r_values = roche_r.(Ω0, lagrange1_x, mass_quotient, vertices(spherical_mesh))
    return GeoTable(
        spherical_mesh,
        Dict(0 => (r = r_values,))
    )
end

function simplexify_quadrangles(mesh::SimpleMesh)
    
end

end