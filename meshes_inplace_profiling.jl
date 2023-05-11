begin
	using Revise
	includet("./Roche_inplace.jl")
	using ..Roche_inplace
end

using Meshes

using Turing
using Profile
using PProf

using BenchmarkTools



begin
	mass_quotients = 0.1:0.3:10

	spherical_mesh = 
		discretize(Sphere((0. ,0., 0.), 1.), RegularDiscretization(16)) |>
		Rotate(Vec(0., 0., 1.), Vec(1., 0., 0.)) |>
		simplexify

	interpolated_roche = InterpolatedRoche(spherical_mesh, mass_quotients, (:f,))
end


m = interpolated_roche(0.5)


f(r) = r / r

apply_radial_function!(m, f, :f)


average_over_faces!(m, :f)


@code_warntype integrate_data_over_triangular_mesh(m, :f, Vec(0., 1., 0.))


@profview for _ in 1:5000
	integrate_data_over_triangular_mesh(m, :f, Vec(0., 1., 0.))
end

@profview for _ in 1:5000
	integrate_data_over_triangular_mesh2(m, :f, Vec(0., 1., 0.))
end



@profview for q ∈ 0.1:0.3:10
	m = interpolated_roche(q)
	apply_radial_function!(m, f, :f)
	integrate_data_over_triangular_mesh2(m, :f, (0., 0., 1.))
end





@model function model(interp_mesh, integral_value)
	q ~ Uniform(0.1, 1.)
	ϕ ~ Uniform(0., π)
	m = interp_mesh(q)
	apply_radial_function!(m, f, :f)
	average_over_faces!(m, :f)
	i = integrate_data_over_triangular_mesh(m, :f, (cos(ϕ), sin(ϕ), 0.))
	integral_value ~ Normal(i, 1e-3)
end


Profile.Allocs.clear()


#Profile.Allocs.@profile sample(model(interpolated_roche, 0.61), NUTS(), 1000)
Profile.Allocs.@profile for q ∈ 0.1:0.3:10, ϕ ∈ 0.:0.1:π
	m = interpolated_roche(q)
	apply_radial_function!(m, f, :f)
	average_over_faces!(m, :f)
	integrate_data_over_triangular_mesh(m, :f, (cos(ϕ), sin(ϕ), 0.))
end

results = Profile.Allocs.fetch()

PProf.Allocs.pprof(results)


@profview for q ∈ 0.1:0.3:10, ϕ ∈ 0.:0.1:π
	m = interpolated_roche(q)
	apply_radial_function!(m, f, :f)
	average_over_faces!(m, :f)
	integrate_data_over_triangular_mesh(m, :f, (cos(ϕ), sin(ϕ), 0.))
end


@profview sample(model(interpolated_roche, 0.61), NUTS(), 1000)

@btime sample(model(interpolated_roche, 0.61), NUTS(), 1000)

@allocated sample(model(interpolated_roche, 0.61), NUTS(), 1000)


@code_warntype average_over_faces!(m, :f)