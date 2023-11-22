begin
    using Revise
    using TuringStars

    using DelimitedFiles
    using DataFrames

    using Turing

    using LombScargle

    using Profile
end

begin
    points = readdlm("stars/T_CrB_JK.dat")[2:end, :]
    points = map(x -> isa(x, Number) ? x : missing, points)
    points = DataFrame(points, [:day, :J, :J_err, :K, :K_err])
end

estimated_period = 227.3200498131982

interpolated_mesh = InterpolatedRocheMesh(64, 0.1:0.1:10)

initial_params = (;
	mass_quotient = 0.5,
	initial_phase = 0.77,
	observer_angle = π/2,
	temperature_at_bottom = 3500,
	offset = 17.25 # 17.17
)

model_params = ModelParams(
	period = estimated_period,
	β = 0.25,
	σ = 0.1,
	luminocity_function = black_body_K,
	temperature_at_bottom = Normal(initial_params.temperature_at_bottom, 500.),
	darkening_function = claret_darkening,
	darkening_coefficients = (1.3113, -1.2998, 1.0144, -0.3272),
	measurements_t = points.day,
	measurements_y = points.K
)

chain_params = ChainParams(
	model_params = model_params,
	n_samples = 10,
	init_params = initial_params,
	sampler = NUTS()
)


@profview samples = TuringStars.ChainCache._uncached_sample(chain_params)