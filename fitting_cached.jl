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
	points = map(x -> isa(x, Number) ? x : NaN, points)
	points = DataFrame(points, [:day, :J, :J_err, :K, :K_err])
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

# ╔═╡ 960ab30d-a1fa-4803-a4d4-d0860286ba87
initial_params = (;
	mass_quotient = 0.5,
	initial_phase = 0.77,
	observer_angle = π/2,
	temperature_at_bottom = 5000,
	offset = 18.14 # 41.4, 18.14, 18.1
)

# ╔═╡ 97fd2129-d706-480c-a97d-9804027d8b40
model_params = ModelParams(
	period = estimated_period,
	β = 0.25,
	fixed_σ = 0.1,
	luminocity_function = black_body_K,
	measurements_t = points.day,
	measurements_y = points.K
)

# ╔═╡ 00044db4-b168-44be-9d39-87d27b7d330d
begin
	scatter(
		points.day .% (estimated_period),
		points.K,
		yerr = points.K_err,
		markersize = 2,
		xlabel = "Julian day % period",
		ylabel = "Звёздная величина",
		title = "K",
		yflip = true
	)

	local days = 0:estimated_period
	phases = @. initial_params.initial_phase + days / estimated_period * 2π
	
	local vals = star_magnitude(
		phases;
		initial_params[(:mass_quotient, :observer_angle, :temperature_at_bottom)]...,
		interpolated_mesh,
		β = model_params.β,
		luminocity_function = model_params.luminocity_function
	)

	vals .+= initial_params.offset

	plot!(days, vals)
end

# ╔═╡ c88314a3-cd9e-42b2-acee-4d613b1b36e1
chain_params = ChainParams(
	model_params = model_params,
	n_samples = 10,
	init_params = initial_params,
	sampler = NUTS()
)

# ╔═╡ eda9134f-b918-42f0-bcfc-e0d601eeeaad
samples = cached_sample(chain_params)

# ╔═╡ a5070b94-48c2-4405-af78-fddd5784161e
chain_params |> JSON3.write |> sha1 |> bytes2hex

# ╔═╡ 174cd8b8-1d1c-4141-a170-6f978f5195e1
begin
	scatter(
		points.day .% estimated_period,
		points.K,
		yerr = points.K_err,
		markersize = 2,
		xlabel = "Julian day % period",
		ylabel = "Звёздная величина",
		title = "K",
		yflip = true
	)

	local days = 0 : estimated_period
	local phases = @. initial_params.initial_phase + days / estimated_period * 2π

	samples_ = sample(samples, 15)
	
	for i in 1 : length(samples_)
	
		local vals = star_magnitude(
			phases;
			mass_quotient = samples_[i][:mass_quotient].data[1],
			observer_angle = samples_[i][:observer_angle].data[1],
			temperature_at_bottom = @something(
				model_params.fixed_temperature_at_bottom,
				samples_[i][:temperature_at_bottom].data[1]
			),
			β = model_params.β,
			interpolated_mesh,
			luminocity_function = model_params.luminocity_function
		)

		vals .+= samples_[i][:offset].data[1]
		plot!(days, vals)
	end
	plot!()
end

# ╔═╡ 45422b39-64d5-4a75-b8c0-8ba0011ba089
plot(samples, bottom_margin = 50Plots.px)

# ╔═╡ Cell order:
# ╠═36db1462-6dbf-11ee-38c4-05d52e2c894c
# ╠═b8bda58e-9ed2-4da0-a67a-6d5990e7389d
# ╠═e28e8c98-caa0-41c0-bb15-53c6679dda6d
# ╠═55c8d8ef-9d4b-4b9c-8838-b91f1f53f8b0
# ╠═5b2930bc-2de0-4388-824a-190d1169cbfe
# ╠═2fe448f3-1744-4bbb-83e7-290a9214e7c8
# ╠═960ab30d-a1fa-4803-a4d4-d0860286ba87
# ╠═00044db4-b168-44be-9d39-87d27b7d330d
# ╠═97fd2129-d706-480c-a97d-9804027d8b40
# ╠═c88314a3-cd9e-42b2-acee-4d613b1b36e1
# ╠═eda9134f-b918-42f0-bcfc-e0d601eeeaad
# ╠═33b862f3-dc6a-46fe-b73e-a7df7af22e92
# ╠═a5070b94-48c2-4405-af78-fddd5784161e
# ╠═174cd8b8-1d1c-4141-a170-6f978f5195e1
# ╠═45422b39-64d5-4a75-b8c0-8ba0011ba089
