module Roche

using Meshes
using Roots
using Interpolations
using LinearAlgebra

export 
    Ω_potential,
    LagrangePoint_X, 
    Ω_critical,
    roche_r,
    StretchToRocheLobe,
    make_roche_mesh,
    make_roche_meshdata,
    InterpolatedRocheMesh,
    integrate_data_over_mesh,
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

    return InterpolatedRocheMesh(spherical_mesh, mass_quotient_knots, interpolants)
end




function (interpolated_mesh::InterpolatedRocheMesh)(mass_quotient)
    points = vertices(interpolated_mesh.spherical_mesh)

    r_list = [
        interpolant(mass_quotient)
        for interpolant ∈ interpolated_mesh.interpolants
    ]
    new_points = [
        Point(coordinates(point) .* r)
        for (point, r) ∈ zip(points, r_list)
    ]
    return meshdata(
        new_points,
        topology(interpolated_mesh.spherical_mesh),
        Dict(:vertices => (r = r_list,))
    )
end



# Functions related to integration


function integrate_data_over_triangular_mesh(mesh_data::MeshData, field_name, direction)
    vertices_values = getfield(values(mesh_data, :vertices), field_name)

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


function apply_radial_function(mesh_data::MeshData, f, field_name)
    r_list = values(mesh_data, :vertices).r
    f_list = f.(r_list)

    return meshdata(
        vertices(domain(mesh_data)),
        topology(domain(mesh_data)),
        Dict(:vertices => (r = r_list, field_name => f_list,))
    ) 
end







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

function make_roche_meshdata(mass_quotient, discretization_method = RegularDiscretization(10))
    mesh = make_roche_mesh(mass_quotient, discretization_method)
    r_list = norm.(coordinates.(vertices(mesh)))
    return meshdata(
        vertices(mesh),
        topology(mesh),
        Dict(:vertices => (r = r_list,))
    )
end

end