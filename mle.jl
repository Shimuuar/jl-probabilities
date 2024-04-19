### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ 11decfb6-8f0a-11ee-1bf8-d3faf8756b8b
begin
	import Pkg
	Pkg.activate(".")

	using Revise
	using TuringStars

	using DelimitedFiles
	using DataFrames

	using Turing
	using Optim

	using Plots
	using StatsPlots
	plotlyjs()
	theme(:juno)

	using LombScargle
end

# ╔═╡ 3d70079b-add5-4b60-9cb2-af7ed5caad87
using Dates

# ╔═╡ 0a209f5e-bd22-4072-9247-c0ddcd42fdb8
begin
	points = readdlm("stars/T_CrB_JK.dat")[2:end, :]
	points = map(x -> isa(x, Number) ? x : missing, points)
	points = DataFrame(points, [:day, :J, :J_err, :K, :K_err])
	points = dropmissing(points)
	points.day .-= points.day[1]
	points
end

# ╔═╡ 18694e43-32c0-424f-ade0-63c6fc8af923
pgram = lombscargle(points.day, points.K)

# ╔═╡ 4b30de9a-e872-422c-9b12-2cc67942aba3
plot(
	freq(pgram),
	power(pgram),
	xlabel = "частота (1/день)",
	title = "Lomb-Scargle periodogram, сигнал K"
)

# ╔═╡ c73161c8-0eb3-40af-b527-40a05e10ae24
begin
	estimated_period = 2findmaxperiod(pgram)[1]
	findmaxperiod(pgram)
end

# ╔═╡ a6dccb45-fa0c-4aff-ae2e-23efc08fd8da
interpolated_mesh = InterpolatedRocheMesh(tetra_sphere(4), 0.1:0.1:10)

# ╔═╡ b782aa83-2e9d-449c-a996-2064c890c025
function plot_line!(model_params, sample, p = missing)
	if p === missing
		p = plot(layout = (length(model_params.channels), 1))
	end

	days = 0 : model_params.period
	phases = days ./ model_params.period .* 2π .+ sample[:initial_phase]

	for (i, (channel, subplot)) ∈ enumerate(zip(model_params.channels, p.subplots))
		vals = star_magnitude(
			phases;
			mass_quotient = sample[:mass_quotient],
			observer_angle = sample[:observer_angle],
			temperature_at_bottom = model_params.temperature_at_bottom,
			interpolated_mesh,
			β = model_params.β,
			luminocity_function = channel.luminocity_function,
			darkening_function = channel.darkening_function,
			darkening_coefs_interpolant = channel.darkening_coefs_interpolant
		)
		vals .+= sample[:offset][i]
		plot!(subplot, days, vals)
	end
	plot!()
end

# ╔═╡ ed6c9a39-4a07-4e98-b893-537d13b1a45d
function plot_garbige(model_params, samples)
	p = plot(
		layout = (2, 1),
		title = ["Спектральный канал K" "Спектральный канал J"],
		legend = false,
		xlabel = ["" "Julian day % period"],
		ylabel = "Звездная величина",
		yflip = true,
		size = (600, 600)
	)

	for (channel, subplot) ∈ zip(model_params.channels, p.subplots)
		t = channel.measurements_t .% model_params.period
		scatter!(
			subplot,
			t,
			channel.measurements_y,
			markersize = 2,
			yerr = channel.σ_measured
		)
	end

	for i ∈ 1 : length(samples)
		sample = samples[i]
		if isa(sample, Chains)
			sample = get_params(sample)
		end
		
		plot_line!(model_params, sample, p)
	end
	plot!()
end

# ╔═╡ 0d54214f-5b15-44cc-bce0-96829c4b5470
channels = [
	ChannelParams(
		measurements_t = points.day,
		measurements_y = points.K,
		darkening_function = claret_darkening,
		darkening_coefs_interpolant = K_coefs_interpolant,
		luminocity_function = black_body_K,
		σ_measured = points.K_err,
		σ_common = FlatPos(0.),
	)
	ChannelParams(
		measurements_t = points.day,
		measurements_y = points.J,
		darkening_function = claret_darkening,
		darkening_coefs_interpolant = J_coefs_interpolant,
		luminocity_function = black_body_J,
		σ_measured = points.J_err,
		σ_common = FlatPos(0.),
	)
]

# ╔═╡ 9d29eca2-57c2-4ace-82de-1472c5509cd8
model_params = ModelParams(
	channels = channels,
	period = estimated_period,
	β = 0.08,
)

# ╔═╡ e0419433-72d5-4d12-9a81-838886dd7317
model = first_model(model_params)

# ╔═╡ b1701d19-4cb7-4c35-943e-6a9cd415c4a1
function logpdf_objective(model)
	var_info = DynamicPPL.VarInfo(model)

	return function f(vars)
		var_info_ = DynamicPPL.unflatten(var_info, vars)
		return logjoint(model, var_info_)
	end
end

# ╔═╡ a3b38672-7e73-4701-abc6-30d057f09754
function MAP_point(model, initial_params)
	var_info = DynamicPPL.VarInfo(model)
	syms = DynamicPPL.syms(var_info)
	x0 = vcat(collect(initial_params[syms])...)

	optim_result = maximize(logpdf_objective(model), x0)
	optim_array = Optim.maximizer(optim_result)

	var_info = DynamicPPL.unflatten(var_info, optim_array)
	return DynamicPPL.values_as(var_info, NamedTuple)
end

# ╔═╡ 991316a4-0a9c-433b-847b-4c4a4e92a24e
initial_params = (;
	mass_quotient = 0.5,
	initial_phase = -1.45,
	observer_angle = π/2 - 0.1,
	temperature_at_bottom = 3500.,
	σ_common = [0.1, 0.1],
	offset = [18.84, 21.15],
)

# ╔═╡ f7c074e7-d5b1-4b2c-9ab6-57e85b237a04
MAP_1 = MAP_point(model, initial_params)

# ╔═╡ 454c33e2-bf88-471d-b51c-6f2010433fef
plot_garbige(model_params, [MAP_1])

# ╔═╡ 66651fe4-7115-40e9-935d-95795757da4e
md"### ``\chi^2``"

# ╔═╡ 28df81cf-e27f-40ee-9ce6-5ff202994dbc
function chi2_value(model_params, MAP_params)

	sum(enumerate(model_params.channels)) do (i, channel)

		phases = channel.measurements_t ./ model_params.period .* 2π .+ MAP_params[:initial_phase]

		predicted_magnitudes = MAP_params[:offset][i]
		predicted_magnitudes = star_magnitude(
			phases;
			mass_quotient = MAP_params[:mass_quotient],
			observer_angle = MAP_params[:observer_angle],
			temperature_at_bottom = model_params.temperature_at_bottom,
			interpolated_mesh,
			β = model_params.β,
			luminocity_function = channel.luminocity_function,
			darkening_function = channel.darkening_function,
			darkening_coefs_interpolant = channel.darkening_coefs_interpolant
		)
		predicted_magnitudes .+= MAP_params[:offset][i]

		σ² = @. channel.σ_measured^2 + MAP_params[Symbol("σ_common[$i]")]

		sum(@. (channel.measurements_y - predicted_magnitudes)^2 / σ²)
	end
end

# ╔═╡ be6c4b94-861e-4907-8a53-dbc3a674664e
chi2_value_ = chi2_value(model_params, MAP_1)

# ╔═╡ bce9145d-23f2-4e22-90bc-0c0089502286
Chi2Dist = Chisq(2 * length(points.day) - 5)

# ╔═╡ 363ea8d3-a174-4186-a967-be659aaa9a72
1 - cdf(Chi2Dist, chi2_value_)

# ╔═╡ eefb61ae-924d-4bb7-9b9a-7235b1f4b826
md"### err = ``f(t)``"

# ╔═╡ b85b3eb7-5257-43d9-9261-0736f7102474
function differences(model_params, MAP_params)

	map(enumerate(model_params.channels)) do (i, channel)

		phases = channel.measurements_t ./ model_params.period .* 2π .+ MAP_params[:initial_phase]

		predicted_magnitudes = MAP_params[:offset][i]
		predicted_magnitudes = star_magnitude(
			phases;
			mass_quotient = MAP_params[:mass_quotient],
			observer_angle = MAP_params[:observer_angle],
			temperature_at_bottom = model_params.temperature_at_bottom,
			interpolated_mesh,
			β = model_params.β,
			luminocity_function = channel.luminocity_function,
			darkening_function = channel.darkening_function,
			darkening_coefs_interpolant = channel.darkening_coefs_interpolant
		)
		predicted_magnitudes .+= MAP_params[:offset][i]

		@. channel.measurements_y - predicted_magnitudes
	end
end

# ╔═╡ c5ad427c-c144-4cc2-b531-ff830beaa21e
plot(
	julian2datetime.(2450199.5 .+ points.day),
	differences(model_params, MAP_1),
	label = ["K" "J"]
)

# ╔═╡ Cell order:
# ╠═11decfb6-8f0a-11ee-1bf8-d3faf8756b8b
# ╠═0a209f5e-bd22-4072-9247-c0ddcd42fdb8
# ╠═18694e43-32c0-424f-ade0-63c6fc8af923
# ╠═4b30de9a-e872-422c-9b12-2cc67942aba3
# ╠═c73161c8-0eb3-40af-b527-40a05e10ae24
# ╠═a6dccb45-fa0c-4aff-ae2e-23efc08fd8da
# ╠═b782aa83-2e9d-449c-a996-2064c890c025
# ╠═ed6c9a39-4a07-4e98-b893-537d13b1a45d
# ╠═0d54214f-5b15-44cc-bce0-96829c4b5470
# ╠═9d29eca2-57c2-4ace-82de-1472c5509cd8
# ╠═e0419433-72d5-4d12-9a81-838886dd7317
# ╠═b1701d19-4cb7-4c35-943e-6a9cd415c4a1
# ╠═a3b38672-7e73-4701-abc6-30d057f09754
# ╠═991316a4-0a9c-433b-847b-4c4a4e92a24e
# ╠═f7c074e7-d5b1-4b2c-9ab6-57e85b237a04
# ╠═454c33e2-bf88-471d-b51c-6f2010433fef
# ╟─66651fe4-7115-40e9-935d-95795757da4e
# ╠═28df81cf-e27f-40ee-9ce6-5ff202994dbc
# ╠═be6c4b94-861e-4907-8a53-dbc3a674664e
# ╠═bce9145d-23f2-4e22-90bc-0c0089502286
# ╠═363ea8d3-a174-4186-a967-be659aaa9a72
# ╟─eefb61ae-924d-4bb7-9b9a-7235b1f4b826
# ╠═b85b3eb7-5257-43d9-9261-0736f7102474
# ╠═c5ad427c-c144-4cc2-b531-ff830beaa21e
# ╠═3d70079b-add5-4b60-9cb2-af7ed5caad87
