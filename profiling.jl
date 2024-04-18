begin
    using Revise
    using TuringStars
	using Meshes

    using DelimitedFiles
    using DataFrames

    using Turing

    using Profile
end

begin
	points = readdlm("stars/T_CrB_JK.dat")[2:end, :]
	points = map(x -> isa(x, Number) ? x : missing, points)
	points = DataFrame(points, [:day, :J, :J_err, :K, :K_err])
	points = dropmissing(points)
	points.day .-= points.day[1]
	points
end

estimated_period = 227.3200498131982

interpolated_mesh = InterpolatedRocheMesh(tetra_sphere(4), 0.1:0.1:10)

channels = [
	ChannelParams(
		measurements_t = points.day,
		measurements_y = points.K,
		darkening_function = claret_darkening,
		darkening_coefs_interpolant = T -> (1.3113, -1.2998, 1.0144, -0.3272),
		luminocity_function = black_body_K,
		σ_measured = points.K_err,
		σ_common = FlatPos(0.),
	)
	ChannelParams(
		measurements_t = points.day,
		measurements_y = points.J,
		darkening_function = claret_darkening,
		darkening_coefs_interpolant = T -> (1.2834, -1.4623, 1.5046, -0.5507),
		luminocity_function = black_body_J,
		σ_measured = points.J_err,
		σ_common = FlatPos(0.),
	)
]

initial_params = (;
	mass_quotient = 0.5,
	initial_phase = -1.45,
	observer_angle = π/2 - 0.1,
	temperature_at_bottom = 3500.,
	σ_common = [0.1, 0.1],
	offset = [18.84, 21.14],
)

model_params = ModelParams(
	channels = channels,
	period = estimated_period,
	β = 0.08,
)

chain_params = ChainParams(
	model_params = model_params,
	n_samples = 20,
	init_params = initial_params,
	sampler = NUTS()
)


@profview @time samples = TuringStars.ChainCache._uncached_sample(chain_params)



mesh = interpolated_mesh(1.0)
mesh = avg_over_faces(mesh, :g)
mesh = apply_function(mesh, g -> abs(g)^0.08, :g, :T)
mesh = apply_function(mesh, T -> T^4, :T, :L)
mesh = apply_function(mesh, T -> (1.2834, -1.4623, 1.5046, -0.5507), :T, :darkening_coefs)

normals = calc_function_on_faces(mesh, normalized_normal)
areas = calc_function_on_faces(mesh, area)

d = (1/√2, 1/√2, 0.)


@profview for _ in 1:1000
	integrate_data_over_mesh(mesh, :g, d, normals, areas, claret_darkening)
end



@profview for _ in 1:100
	star_magnitude(
		0.0:0.01:2π,
		mass_quotient=0.5,
		observer_angle=π/2 - 0.1,
		temperature_at_bottom=3500.,
		β=0.08,
		interpolated_mesh=interpolated_mesh,
		luminocity_function=black_body_K,
		darkening_function=claret_darkening,
		darkening_coefs_interpolant = K_coefs_interpolant
	)
end
