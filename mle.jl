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
	theme(:default)

	using LombScargle
end

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
interpolated_mesh = InterpolatedRocheMesh(64, 0.1:0.1:10)

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
			darkening_coefficients = channel.darkening_coefficients
		)
		vals .+= sample[Symbol("offset[$i]")]
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
		darkening_coefficients = (1.3113, -1.2998, 1.0144, -0.3272),
		luminocity_function = black_body_K,
		σ_measured = points.K_err,
		σ_common = FlatPos(0.),
	)
	ChannelParams(
		measurements_t = points.day,
		measurements_y = points.J,
		darkening_function = claret_darkening,
		darkening_coefficients = (1.2834, -1.4623, 1.5046, -0.5507),
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

# ╔═╡ 218df127-0f16-41b2-89d4-02be6b331936
mle_estimate = optimize(model, MLE())

# ╔═╡ 2356ab8e-2098-48eb-abfa-7371f2b7a6ee
plot_garbige(model_params, [mle_estimate.values])

# ╔═╡ a071af8e-50cc-40f2-96dd-49c13266a885
mle_estimate.values[:mass_quotient]

# ╔═╡ 90bb0d33-548e-4c33-ba8c-b11926529a9c
model_params2 = ModelParams(
	channels = channels,
	period = estimated_period,
	β = 0.25,
)

# ╔═╡ 995f20d1-90e1-4d74-a7e2-99c744e8eb24
model2 = first_model(model_params2)

# ╔═╡ d1cdc25b-1d80-41a3-8561-9c8e187ebca3
mle_estimate2 = optimize(model2, MLE())

# ╔═╡ b1c4696f-e4c6-4006-aab6-0745879f2131
mle_estimate2.values[:initial_phase]

# ╔═╡ 53894ab6-e2c7-49cb-b6f6-7159f9bfd306
mle_estimate2.values[:mass_quotient]

# ╔═╡ 1bfa3efb-9a5c-45bd-9a79-6e2093b57f1d
plot_garbige(model_params2, [mle_estimate2.values])

# ╔═╡ 45afe5bd-fcfe-4d67-bd14-d9d607bf8724
plot(points.day ./ 365)

# ╔═╡ f5abbcbf-4dc7-46ba-bf98-7ad0dd30cb99
function plot_rubbish2(model_params, sample)
	p = plot(
		layout = (2, 1),
		title = ["Спектральный канал K" "Спектральный канал J"],
		legend = false,
		xlabel = ["" "фаза"],
		ylabel = "Звездная величина",
		yflip = true,
		size = (600, 600),
		xlim = (-0.1, 1.1),
		margin = 12Plots.px,
	)

	initial_phase = sample[:initial_phase] / 2π

	for (channel, subplot) ∈ zip(model_params.channels, p.subplots)
		phases = (channel.measurements_t .% model_params.period) ./ model_params.period .+ initial_phase
		scatter!(
			subplot,
			phases,
			channel.measurements_y,
			markersize = 2,
			yerr = channel.σ_measured,
			markercolor = 1,
		)
		scatter!(
			subplot,
			phases .+ 1,
			channel.measurements_y,
			markersize = 2,
			yerr = channel.σ_measured,
			markercolor = 1,
		)
		scatter!(
			subplot,
			phases .-1,
			channel.measurements_y,
			markersize = 2,
			yerr = channel.σ_measured,
			markercolor = 1,
		)
	end

	phases = -0.1 : 0.01 : 1.1

	for (c, (channel, subplot)) ∈ enumerate(zip(model_params.channels, p.subplots))

		if isa(sample, Chains)
			sample = get_params(sample)
		end

		vals = star_magnitude(
			phases .* 2π;
			mass_quotient = sample[:mass_quotient],
			observer_angle = sample[:observer_angle],
			temperature_at_bottom = model_params.temperature_at_bottom,
			interpolated_mesh,
			β = model_params.β,
			luminocity_function = channel.luminocity_function,
			darkening_function = channel.darkening_function,
			darkening_coefficients = channel.darkening_coefficients
		)
		vals .+= sample[Symbol("offset[$c]")]

		plot!(
			subplot,
			phases,
			vals,
			color = 2
		)
	end
	p
end

# ╔═╡ 6373adaa-9a74-42a2-b159-f9190872bc83
mle_estimate.values[Symbol("offset[1]")]

# ╔═╡ 42837d09-13a9-46df-8f65-fb8f241f6fba
plot_rubbish2(model_params, mle_estimate.values)

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
# ╠═218df127-0f16-41b2-89d4-02be6b331936
# ╠═2356ab8e-2098-48eb-abfa-7371f2b7a6ee
# ╠═a071af8e-50cc-40f2-96dd-49c13266a885
# ╠═90bb0d33-548e-4c33-ba8c-b11926529a9c
# ╠═995f20d1-90e1-4d74-a7e2-99c744e8eb24
# ╠═d1cdc25b-1d80-41a3-8561-9c8e187ebca3
# ╠═b1c4696f-e4c6-4006-aab6-0745879f2131
# ╠═53894ab6-e2c7-49cb-b6f6-7159f9bfd306
# ╠═1bfa3efb-9a5c-45bd-9a79-6e2093b57f1d
# ╠═45afe5bd-fcfe-4d67-bd14-d9d607bf8724
# ╠═f5abbcbf-4dc7-46ba-bf98-7ad0dd30cb99
# ╠═6373adaa-9a74-42a2-b159-f9190872bc83
# ╠═42837d09-13a9-46df-8f65-fb8f241f6fba
