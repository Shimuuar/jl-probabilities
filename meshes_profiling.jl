begin
    using Revise
    includet("./Roche.jl")
    using ..Roche
end

using Meshes

using Turing
using Profile

using BenchmarkTools



begin
    mass_quotients = 0.1:0.3:10

    spherical_mesh =
        discretize(Sphere((0. ,0., 0.), 1.), RegularDiscretization(16)) |>
        Rotate(Vec(0., 0., 1.), Vec(1., 0., 0.)) |>
        simplexify

    interpolated_mesh = InterpolatedRocheMesh(spherical_mesh, mass_quotients)
end


m = interpolated_mesh(0.5)


f(r) = 1.

m2 = apply_radial_function(m, f, :f)



@code_warntype integrate_data_over_triangular_mesh(m2, :f, Vec(0., 1., 0.))


@profview for _ in 1:5000
    integrate_data_over_triangular_mesh(m2, :f, Vec(0., 1., 0.))
end




@model function model(interp_mesh, integral_value)
    q ~ Uniform(0.1, 1.)
    ϕ ~ Uniform(0., π)
    m = interp_mesh(q)
    m = apply_radial_function(m, f, :f)
    i = integrate_data_over_triangular_mesh(m, :f, (cos(ϕ), sin(ϕ), 0.))
    integral_value ~ Normal(i, 1e-3)
end


@profview sample(model(interpolated_mesh, 0.61), NUTS(), 1000)

@btime sample(model(interpolated_mesh, 0.61), NUTS(), 1000)

