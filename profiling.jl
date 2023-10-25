begin
    using Revise
    using TuringStars

    using DelimitedFiles
    using DataFrames

    using Turing

    using Plots
    using StatsPlots
    plotlyjs()
    theme(:juno)

    using LombScargle

    using Profile
end

begin
    points = readdlm("stars/T_CrB_JK.dat")[2:end, :]
    points = map(x -> isa(x, Number) ? x : NaN, points)
    points = DataFrame(points, [:day, :J, :J_err, :K, :K_err])
end

estimated_period = 227.3200498131982

interpolated_mesh = InterpolatedRocheMesh(64, 0.1:0.1:10)

initial_params = (;
	mass_quotient = 0.5,
	initial_phase = 0.77,
	observer_angle = π/2,
	temperature_at_bottom = 5000,
	offset = 18.14 # 41.4, 18.14, 18.1
)

model_params = ModelParams(
	period = estimated_period,
	β = 0.25,
	fixed_σ = 0.1,
	luminocity_function = black_body_K,
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