### A Pluto.jl notebook ###
# v0.19.29

using Markdown
using InteractiveUtils

# ╔═╡ 36db1462-6dbf-11ee-38c4-05d52e2c894c
begin
	import Pkg
	Pkg.activate(".")

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
end

# ╔═╡ 33b862f3-dc6a-46fe-b73e-a7df7af22e92
using JSON3, SHA

# ╔═╡ b8bda58e-9ed2-4da0-a67a-6d5990e7389d
begin
	points = readdlm("stars/T_CrB_JK.dat")[2:end, :]
	points = map(x -> isa(x, Number) ? x : missing, points)
	points = DataFrame(points, [:day, :J, :J_err, :K, :K_err])
	points = dropmissing(points)
	points.day .-= points.day[1]
	points
end

# ╔═╡ e28e8c98-caa0-41c0-bb15-53c6679dda6d
pgram = lombscargle(points.day, points.K)

# ╔═╡ 55c8d8ef-9d4b-4b9c-8838-b91f1f53f8b0
plot(
	freq(pgram),
	power(pgram),
	xlabel = "частота (1/день)",
	title = "Lomb-Scargle periodogram, сигнал K"
)

# ╔═╡ 5b2930bc-2de0-4388-824a-190d1169cbfe
begin
	estimated_period = 2findmaxperiod(pgram)[1]
	findmaxperiod(pgram)
end

# ╔═╡ 2fe448f3-1744-4bbb-83e7-290a9214e7c8
interpolated_mesh = InterpolatedRocheMesh(64, 0.1:0.1:10)

# ╔═╡ d9b9b851-cca0-4a40-a627-7dec9c5da6c1
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
		vals .+= sample[:offset][i]
		plot!(subplot, days, vals)
	end
	plot!()
end

# ╔═╡ 275cb92f-d5d1-4fb9-acb3-5d2317f84a2b
function plot_garbige(model_params, samples)
	p = plot(
		layout = (2, 1),
		title = ["K" "J"],
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

# ╔═╡ 960ab30d-a1fa-4803-a4d4-d0860286ba87
initial_params = (;
	mass_quotient = 0.5,
	initial_phase = -1.45,
	observer_angle = π/2 - 0.1,
	temperature_at_bottom = 3500.,
	σ_common = [0.1, 0.1],
	offset = [18.84, 21.15],
)

# ╔═╡ 30a74a85-c431-469c-bf3d-00190db36c56
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
		σ_measured = points.K_err,
		σ_common = FlatPos(0.),
	)
]

# ╔═╡ 97fd2129-d706-480c-a97d-9804027d8b40
model_params = ModelParams(
	channels = channels,
	period = estimated_period,
	β = 0.08,
)

# ╔═╡ 00044db4-b168-44be-9d39-87d27b7d330d
plot_garbige(model_params, [initial_params])

# ╔═╡ c88314a3-cd9e-42b2-acee-4d613b1b36e1
chain_params = ChainParams(
	model_params = model_params,
	n_samples = 2000,
	init_params = initial_params,
	sampler = NUTS()
)

# ╔═╡ eda9134f-b918-42f0-bcfc-e0d601eeeaad
samples = cached_sample(chain_params)

# ╔═╡ a5070b94-48c2-4405-af78-fddd5784161e
chain_params |> JSON3.write |> sha1 |> bytes2hex

# ╔═╡ 94abc73f-8f2e-42f5-86d2-17836d645ec2
macro get_number(symbol, parameters, sample)
	@eval isa($parameters.$symbol, Number) ? $parameters.$symbol : $sample[symbol].data[1]
end

# ╔═╡ 174cd8b8-1d1c-4141-a170-6f978f5195e1
begin
	samples_ = sample(samples, 15)
	plot_garbige(model_params, samples_)
end

# ╔═╡ 45422b39-64d5-4a75-b8c0-8ba0011ba089
plot(samples, bottom_margin = 50Plots.px)

# ╔═╡ 61b8f08b-4252-4ea7-9b27-37771331de77


# ╔═╡ Cell order:
# ╠═36db1462-6dbf-11ee-38c4-05d52e2c894c
# ╠═b8bda58e-9ed2-4da0-a67a-6d5990e7389d
# ╠═e28e8c98-caa0-41c0-bb15-53c6679dda6d
# ╠═55c8d8ef-9d4b-4b9c-8838-b91f1f53f8b0
# ╠═5b2930bc-2de0-4388-824a-190d1169cbfe
# ╠═2fe448f3-1744-4bbb-83e7-290a9214e7c8
# ╠═d9b9b851-cca0-4a40-a627-7dec9c5da6c1
# ╠═275cb92f-d5d1-4fb9-acb3-5d2317f84a2b
# ╠═960ab30d-a1fa-4803-a4d4-d0860286ba87
# ╠═00044db4-b168-44be-9d39-87d27b7d330d
# ╠═30a74a85-c431-469c-bf3d-00190db36c56
# ╠═97fd2129-d706-480c-a97d-9804027d8b40
# ╠═c88314a3-cd9e-42b2-acee-4d613b1b36e1
# ╠═eda9134f-b918-42f0-bcfc-e0d601eeeaad
# ╠═33b862f3-dc6a-46fe-b73e-a7df7af22e92
# ╠═a5070b94-48c2-4405-af78-fddd5784161e
# ╠═94abc73f-8f2e-42f5-86d2-17836d645ec2
# ╠═174cd8b8-1d1c-4141-a170-6f978f5195e1
# ╠═45422b39-64d5-4a75-b8c0-8ba0011ba089
# ╠═61b8f08b-4252-4ea7-9b27-37771331de77
