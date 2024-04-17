module Roche

using Meshes
using GeoTables
using Roots
using Interpolations
using LinearAlgebra
using ForwardDiff

export
    Ω_potential,
    Ω_grad,
    LagrangePoint_X,
    Ω_critical,
    roche_r,
    InterpolatedRocheMesh,
    tetra_sphere,
    integrate_data_over_mesh,
    apply_function,
    avg_over_faces,
    calc_function_on_faces,
    luminocity_at_point,
    normalized_normal


"""
Формула 6.1 Черепащука (стр 88)
"""
function Ω_potential(r::Number; mass_quotient::Number, point_on_unit_sphere::Point)
    λ, μ, ν = coordinates(point_on_unit_sphere)
    return 1/r +
        mass_quotient * (1 / √(1 + r^2 - 2r*λ) - λ*r) +
        (1. + mass_quotient) / 2 * r^2 * (1. - ν^2)
end



"""
Переделанная формула 6.1 Черепащука, выраженная через x,y,z вместо r,λ,μ,ν

Нужна для вычисления градиента
"""
function Ω_potential(mass_quotient, coords)
    x, y, z = coords
    r = hypot(x, y, z)
    return 1/r +
        mass_quotient * (1 / √(1 + r^2 - 2x) - x) +
        (1. + mass_quotient) / 2 * (x^2 + y^2)
end

function Ω_grad(mass_quotient, coords::Vec)
    return ForwardDiff.gradient(
        xyz -> Ω_potential(mass_quotient, xyz),
        coords
    )
end



"""
Ищем точку Лагранжа по оси X, соединяющей две звезды.

Для этого приравниваем ∂Ω/∂x к нулю и ищем корень на отрезке (0, 1).

Формула для производной — 6.2 Черепащука (стр 88)
"""
function LagrangePoint_X(mass_quotient)
    function Ω_x_derivative(x)
        -1 / x^2 - 
        mass_quotient * (1. - 1 / (x-1)^2) + 
        (1 + mass_quotient) * x
    end
    return find_zero(Ω_x_derivative, (0., 1.))
end



"""
Значение Ω-потенциала в точке Лагранжа (и вообще в любой точке на поверхности полости Роша)
"""
function Ω_critical(lagrange1_x, mass_quotient)
    return Ω_potential(lagrange1_x; mass_quotient, point_on_unit_sphere=Point(1., 0., 0.))
end



"""
Радиус полости Роша вдоль направления `point_on_unit_sphere`.

Численно решаем уравнение f(r) = Ω(r, λ, μ, ν) - Ω0 = 0, где Ω0 — значение Ω-потенциала в точке Лагранжа.
"""
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

"""
Интерполированная по соотношению масс сетка полости Роша.

Для InterpolatedRocheMesh переопределен оператор вызова функции. Если применить её к соотношению масс, то получится Meshes.GeoTable, содержащая сетку полости Роша при данном соотношении масс и значения |g| в каждой точке, нормированные на значение |g| в точке (x, y, z) = (-1, 0, 0).
"""
struct InterpolatedRocheMesh
    spherical_mesh
    mass_quotient_knots
    r_interpolants
    g_interpolants
end



"""
Конструктор InterpolatedRocheMesh, если уже построена сферическая сетка.
"""
function InterpolatedRocheMesh(spherical_mesh::SimpleMesh, mass_quotient_knots)

    r_values = zeros(nvertices(spherical_mesh), length(mass_quotient_knots))
    g_values = zeros(nvertices(spherical_mesh), length(mass_quotient_knots))

    for (mass_quotient, r_values_for_quotient, g_values_for_quotient) ∈
            zip(mass_quotient_knots, eachcol(r_values), eachcol(g_values))

        lagrange1_x = LagrangePoint_X(mass_quotient)
        Ω0 = Ω_critical(lagrange1_x, mass_quotient)
        r_values_for_quotient .= roche_r.(Ω0, lagrange1_x, mass_quotient, vertices(spherical_mesh))

        p0 = Point(-1., 0., 0.)
        r0 = roche_r(Ω0, lagrange1_x, mass_quotient, p0)
        g0 = norm(Ω_grad(mass_quotient, coordinates(p0) .* r0))

        coords = coordinates.(vertices(spherical_mesh)) .* r_values_for_quotient
        g_values_for_quotient .= norm.(Ω_grad.(mass_quotient, coords)) ./ g0
        r_values_for_quotient ./= r0
    end

    r_interpolants = [
        cubic_spline_interpolation(mass_quotient_knots, r_values_for_vertex)
        for r_values_for_vertex ∈ eachrow(r_values)
    ]

    g_interpolants = [
        cubic_spline_interpolation(mass_quotient_knots, g_values_for_vertex)
        for g_values_for_vertex ∈ eachrow(g_values)
    ]

    return InterpolatedRocheMesh(spherical_mesh, mass_quotient_knots, r_interpolants, g_interpolants)
end


"""
Создаёт сферическую сетку, используя алгоритм Катмулла-Кларка.
Начальная сетка — тетраэдр.
"""
function tetra_sphere(catmullclark_iterations)
    box = Tetrahedron(
        (-√2/3, -√2/√3, -1/3),
        (0, 0, 1),
        (2√2/3, 0, -1/3),
        (-√2/3, √2/√3, -1/3),
    ) |> boundary |> discretize

    for _ ∈ 1 : catmullclark_iterations
        box = refine(box, CatmullClark())
    end

    points = map(vertices(box)) do point
        coords = coordinates(point)
        coords = coords ./ norm(coords)
        Point(coords...)
    end

    return SimpleMesh(points, topology(box)) |> simplexify
end


"""
Meshes.GeoTable, содержащая:
* сетку полости Роша при данном соотношении масс,
* значения |g| в каждой точке, нормированные на значение |g| в точке (x, y, z) = (-1, 0, 0),
* значения T = abs(g)^β в каждой грани,

"""
function (interpolated_mesh::InterpolatedRocheMesh)(mass_quotient)
    points = vertices(interpolated_mesh.spherical_mesh)

    r_list = mass_quotient .|> interpolated_mesh.r_interpolants
    g_list = mass_quotient .|> interpolated_mesh.g_interpolants

    new_points = [
        Point(coordinates(point) .* r...)
        for (point, r) ∈ zip(points, r_list)
    ]
    new_mesh = SimpleMesh(
        new_points,
        topology(interpolated_mesh.spherical_mesh)
    )
    return GeoTable(
        new_mesh,
        Dict(0 => (g = g_list,))
    )
end



# Functions related to integration

"""
Интеграл по сетке.
Аргументы:
- geo_table::Meshes.GeoTable — сетка; должна хранить коэффициенты потемнения к краю для каждой грани: 2 => :darkening_coefs
- field_name — имя поля, под которым в гранях geo_table хранятся значения, которые нужно интегрировать
- direction — направление на наблюдателя
- normals — предварительно вычисленные внешние нормали к граням сетки
- areas — предварительно вычисленные площади граней
- darkening_function — функция потемнения к краю, вызываемая как darkening_function(cosine, darkening_coefficients...)
"""
function integrate_data_over_mesh(geo_table::GeoTable, field_name, direction, normals, areas,
                                  darkening_function)
    faces_values = getfield(values(geo_table, 2), field_name)
    darkening_coefs = getfield(values(geo_table, 2), :darkening_coefs)

    sum(zip(faces_values, darkening_coefs, normals, areas)) do (val, darkening_coef, n, A)
        cosine = n ⋅ direction
        if cosine < 0
            return zero(cosine)
        else
            darkening = darkening_function(cosine, darkening_coef...)
            return cosine * A * val * darkening
        end
    end
end


"""
Применяет функцию к значениям, которые лежат в гранях GeoTable, и возвращает новую GeoTable c добавленным полем результатов.
"""
function apply_function(geo_table::GeoTable, f, arg_field_name, result_field_name, dim=2)
    data = getfield(geo_table, :values)

    f_list = f.(getfield(data[dim], arg_field_name))

    data = merge_at(dim, data, result_field_name, f_list)

    return GeoTable(domain(geo_table), data)
end

"""
Добавить Dict(dim => (field_name = array,)) к data.
"""
function merge_at(dim, data, field_name, array)
    data_dim = get(data, dim, NamedTuple())
    data_dim = merge(
        data_dim,
        NamedTuple{(field_name,)}((array,))
    )
    data = merge(data, Dict(dim => data_dim))
    return data
end


"""
Вспомогательная функция, чтобы усреднить значения в вершинах грани.
Аргументы:
- vertices_values — массив значений во всех вершинах
- connection::Meshes.Connectivity — индексы нескольких вершин, образующих грань
"""
function avg_over_face(vertices_values, connection)
    index = indices(connection)
    return sum(vertices_values[i] for i ∈ index) / length(index)
end


"""
Усредняет значения по граням и записывает их в GeoTable.
"""
function avg_over_faces(geo_table::GeoTable, field_name, dim=2)
    vertices_values = getfield(values(geo_table, 0), field_name)
    connections = faces(topology(domain(geo_table)), dim)
    faces_values = map(connections) do connection
        avg_over_face(vertices_values, connection)
    end
    data = getfield(geo_table, :values)
    data = merge_at(dim, data, field_name, faces_values)
    return GeoTable(domain(geo_table), data)
end

"""
Применяет функцию к граням и возвращает массив значений. Грани имеют тип, например, Meshes.Triangle{3, Float64}.
"""
function calc_function_on_faces(geo_table::GeoTable, f, dim=2)
    f.(faces(domain(geo_table), dim))
end


"""
Единичная нормаль к грани.
Эта функция нужна начиная с версии Meshes.jl 0.39.0, в которой функция `normal` стала возвращать ненормированную нормаль.
"""
function normalized_normal(face::Meshes.Ngon{N, T}) where {N, T}
    n = normal(face)
    return n ./ norm(n)
end

end