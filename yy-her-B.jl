### A Pluto.jl notebook ###
# v0.19.42

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
	#plotlyjs()
	#theme(:juno)

	using LombScargle
end

# ╔═╡ 234e80df-af67-44ad-8f06-3c3d403dcd25
using KernelDensity

# ╔═╡ 360119a3-9cba-4e40-ad6c-88a0be627b1e
using Roots

# ╔═╡ 68c78f49-b573-4f29-96df-ca65f5ea5ae9


# ╔═╡ b8bda58e-9ed2-4da0-a67a-6d5990e7389d
begin
	points = readdlm("stars/Yy_her.dat")[:, 2:end]
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
#estimated_period = 586.75

# ╔═╡ 7e8b804f-b511-4359-8c44-286743a2ff9b
estimated_period

# ╔═╡ 2fe448f3-1744-4bbb-83e7-290a9214e7c8
interpolated_mesh = InterpolatedRocheMesh(tetra_sphere(4), 0.1:0.1:10)

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
			darkening_coefs_interpolant = channel.darkening_coefs_interpolant
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
		title = ["Спектральный канал K" "Спектральный канал J"],
		legend = false,
		xlabel = ["" "Julian day % period"],
		ylabel = "Звездная величина",
		yflip = true,
		size = (600, 600),
		margin = 12Plots.px,
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

# ╔═╡ 232ace15-13a9-4afe-9468-d6e54d796470
function plot_rubbish(model_params, samples)
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

	initial_phase = mean(samples[:initial_phase]) / 2π

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
		vals = Array{Float64}(undef, length(samples), length(phases))

		for s ∈ 1 : length(samples)
			sample = samples[s]
			if isa(sample, Chains)
				sample = get_params(sample)
			end

			vals[s, :] = star_magnitude(
				phases .* 2π;
				mass_quotient = sample[:mass_quotient],
				observer_angle = sample[:observer_angle],
				temperature_at_bottom = model_params.temperature_at_bottom,
				interpolated_mesh,
				β = model_params.β,
				luminocity_function = channel.luminocity_function,
				darkening_function = channel.darkening_function,
				darkening_coefs_interpolant = channel.darkening_coefs_interpolant
			)
			vals[s, :] .+= sample[:offset][c]
		end

		means = mean(vals, dims = 1)[1, :]
		stds = std(vals, dims = 1)[1, :]

		plot!(
			subplot,
			phases,
			means,
			alpha = 0,
			ribbon = stds,
			color = 1
		)
	end
	p
end

# ╔═╡ 960ab30d-a1fa-4803-a4d4-d0860286ba87
initial_params = (;
	mass_quotient = 0.5,
	initial_phase = -.7,
	observer_angle = π/2 - 0.1,
	σ_common = [0.1, 0.1],
	offset = [22.0, 24.3],
)

# ╔═╡ 30a74a85-c431-469c-bf3d-00190db36c56
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

# ╔═╡ a5fa116c-46c3-49fb-93ad-6cddf46bea04
plot_garbige(
	ModelParams(
		channels = channels,
		period = estimated_period + 0.5,
		β = 0.08,
	),
	[]
)

# ╔═╡ 97fd2129-d706-480c-a97d-9804027d8b40
model_params = ModelParams(
	channels = channels,
	period = estimated_period,
	β = 0.08,
)

# ╔═╡ e4d38acb-9327-4f40-84d3-a33e56bea394
plot_garbige(model_params, [])

# ╔═╡ 00044db4-b168-44be-9d39-87d27b7d330d
plot_garbige(model_params, [initial_params])

# ╔═╡ c88314a3-cd9e-42b2-acee-4d613b1b36e1
chain_params = ChainParams(
	model_params = model_params,
	n_samples    = 10000,
	n_chains     = 8,
	init_params  = initial_params,
	sampler      = NUTS()
)

# ╔═╡ eda9134f-b918-42f0-bcfc-e0d601eeeaad
samples = cached_sample(chain_params)

# ╔═╡ 174cd8b8-1d1c-4141-a170-6f978f5195e1
begin
	samples_ = sample(samples, 15)
	plot_garbige(model_params, samples_)
end

# ╔═╡ 15ae4b29-3d14-4ad1-811d-de973095f25d
plot_rubbish(model_params, samples)

# ╔═╡ 45422b39-64d5-4a75-b8c0-8ba0011ba089
plot(samples, margin = 10Plots.px, bottom_margin = 50Plots.px)

# ╔═╡ 43973ad5-74c8-4bb7-92c2-372e6fb722dd
density(
	1 ./ samples[:mass_quotient],
	title = "Масса гиганта / масса карлика",
	legend = false,
	size = (600, 300),
	xlim = (0, 2)
)

# ╔═╡ 8ec7dad5-b377-43df-8a7f-ad8c9179d7e2
density(
	samples[:initial_phase] ./ 2π,
	title = "Начальная фаза (в долях периода)",
	legend = false,
	size = (600, 300),
)

# ╔═╡ 9a192e6e-f8ac-48fc-a7c3-7cbc5ebe7554
density(
	samples[:observer_angle] ./ π .* 180,
	title = "Наклонение (°)",
	legend = false,
	size = (600, 300),
)

# ╔═╡ 61b8f08b-4252-4ea7-9b27-37771331de77
scatter(
	1 ./ samples[:mass_quotient],
	samples[:observer_angle] / pi * 180,
	markersize = 0.2,
	xlabel = "Масса гиганта / масса карлика",
	ylabel = "Наклонение (рад.)",
	legend = false
)

# ╔═╡ 52b94bf6-3dc2-46ca-bc70-ad71f47085b2
samples[:observer_angle]

# ╔═╡ fcff208a-6eea-475b-8cee-908457bbc1d3
k = kde((
	1 ./ samples[:mass_quotient].data[:, 1],
	samples[:observer_angle].data[:, 1] .* 180 ./ π
),
	npoints=(256,512),
	#boundary=((0.,8.), (30.0,90.0))
)

# ╔═╡ 46a12905-4167-4f7e-92a6-6a633903f7e9
function get_threshold(kde, confidence_level)
	total = sum(kde.density)
	function percentage(threshold)
		return sum(filter(pdf -> pdf > threshold, kde.density)) / total
	end
	return find_zero(
		threshold -> percentage(threshold) - confidence_level,
		extrema(kde.density),
	)
end

# ╔═╡ dedc4007-cf2e-49ac-860f-b6b27cd2be05
labels_ = [0.997, 0.95, 0.68]

# ╔═╡ 230396d1-c1ce-4b71-a533-ceb745392393
text_lavels = @. string((labels_ * 100)) * "%"

# ╔═╡ f5bccc4e-aefe-4bc7-a6c2-da3de047b5d4
thresholds = get_threshold.(Ref(k), labels_)

# ╔═╡ 4c25f206-eb15-46a5-8685-d299d71695b6
heatmap(k.x, k.y, transpose(k.density))

# ╔═╡ eb82528e-657b-4d80-87c4-5feb9dff95ef
begin
	#theme(:)
	contour(k.x, k.y, transpose(k.density), 
		levels=thresholds,
		xlabel = "Масса гиганта / масса карлика",
		ylabel = "Наклонение",
	)
end

# ╔═╡ 1711c69b-70d1-4163-834b-0d55bf441b87
length(k.x)

# ╔═╡ 58326fc1-8eba-4dba-9071-685129564306
k.y

# ╔═╡ 5d7caa39-2328-4e12-805b-3438ce74faa8
# ╠═╡ disabled = true
#=╠═╡
import PyPlot
  ╠═╡ =#

# ╔═╡ 072d2ba5-4bf8-4d9d-8165-38bcbfe8cd82
# ╠═╡ disabled = true
#=╠═╡
PyPlot.svg(true)
  ╠═╡ =#

# ╔═╡ 58497ce3-5b59-4473-b3c3-9ddea6de1cb9
#=╠═╡
begin
	PyPlot.gcf().clear()

	PyPlot.gca().set_xlabel("Масса гиганта / масса карлика")
	PyPlot.gca().set_ylabel("Наклонение (°)")

	local cs = PyPlot.contour(
		k.x, k.y, k.density',
		levels = thresholds
	)
	PyPlot.clabel(cs, fmt = Dict(zip(
		thresholds,
		text_lavels
	)))
	PyPlot.scatter([max_point[1]], [max_point[2]])
	PyPlot.gcf()
end
  ╠═╡ =#

# ╔═╡ 4693a488-e017-46a8-b077-cffd3b20303e
# ╠═╡ disabled = true
#=╠═╡
function i(q; m = 1.44)
	rad2deg(asin(cbrt(0.322408 / m * (1 + q)^2)))
end
  ╠═╡ =#

# ╔═╡ 3b39b5a4-6cb1-4b80-9d13-730f6a797fd8
# ╠═╡ disabled = true
#=╠═╡
using Optim
  ╠═╡ =#

# ╔═╡ fc855167-775a-4688-adbc-93625d72c3cd
# ╠═╡ disabled = true
#=╠═╡
mx = maximize(x -> pdf(k, x[1], x[2]), [1., 60.])
  ╠═╡ =#

# ╔═╡ a5f25a95-d00f-4cc9-b93e-f02efc812e9d
# ╠═╡ disabled = true
#=╠═╡
max_point = Optim.maximizer(mx)
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═36db1462-6dbf-11ee-38c4-05d52e2c894c
# ╠═68c78f49-b573-4f29-96df-ca65f5ea5ae9
# ╠═b8bda58e-9ed2-4da0-a67a-6d5990e7389d
# ╠═e28e8c98-caa0-41c0-bb15-53c6679dda6d
# ╠═55c8d8ef-9d4b-4b9c-8838-b91f1f53f8b0
# ╠═5b2930bc-2de0-4388-824a-190d1169cbfe
# ╠═7e8b804f-b511-4359-8c44-286743a2ff9b
# ╠═e4d38acb-9327-4f40-84d3-a33e56bea394
# ╠═2fe448f3-1744-4bbb-83e7-290a9214e7c8
# ╠═d9b9b851-cca0-4a40-a627-7dec9c5da6c1
# ╠═275cb92f-d5d1-4fb9-acb3-5d2317f84a2b
# ╠═232ace15-13a9-4afe-9468-d6e54d796470
# ╠═960ab30d-a1fa-4803-a4d4-d0860286ba87
# ╠═00044db4-b168-44be-9d39-87d27b7d330d
# ╠═a5fa116c-46c3-49fb-93ad-6cddf46bea04
# ╠═30a74a85-c431-469c-bf3d-00190db36c56
# ╠═97fd2129-d706-480c-a97d-9804027d8b40
# ╠═c88314a3-cd9e-42b2-acee-4d613b1b36e1
# ╠═eda9134f-b918-42f0-bcfc-e0d601eeeaad
# ╠═174cd8b8-1d1c-4141-a170-6f978f5195e1
# ╠═15ae4b29-3d14-4ad1-811d-de973095f25d
# ╠═45422b39-64d5-4a75-b8c0-8ba0011ba089
# ╠═43973ad5-74c8-4bb7-92c2-372e6fb722dd
# ╠═8ec7dad5-b377-43df-8a7f-ad8c9179d7e2
# ╠═9a192e6e-f8ac-48fc-a7c3-7cbc5ebe7554
# ╠═61b8f08b-4252-4ea7-9b27-37771331de77
# ╠═52b94bf6-3dc2-46ca-bc70-ad71f47085b2
# ╠═234e80df-af67-44ad-8f06-3c3d403dcd25
# ╠═360119a3-9cba-4e40-ad6c-88a0be627b1e
# ╠═fcff208a-6eea-475b-8cee-908457bbc1d3
# ╠═46a12905-4167-4f7e-92a6-6a633903f7e9
# ╠═dedc4007-cf2e-49ac-860f-b6b27cd2be05
# ╠═230396d1-c1ce-4b71-a533-ceb745392393
# ╠═f5bccc4e-aefe-4bc7-a6c2-da3de047b5d4
# ╠═4c25f206-eb15-46a5-8685-d299d71695b6
# ╠═eb82528e-657b-4d80-87c4-5feb9dff95ef
# ╠═1711c69b-70d1-4163-834b-0d55bf441b87
# ╠═58326fc1-8eba-4dba-9071-685129564306
# ╠═5d7caa39-2328-4e12-805b-3438ce74faa8
# ╠═072d2ba5-4bf8-4d9d-8165-38bcbfe8cd82
# ╠═58497ce3-5b59-4473-b3c3-9ddea6de1cb9
# ╠═4693a488-e017-46a8-b077-cffd3b20303e
# ╠═3b39b5a4-6cb1-4b80-9d13-730f6a797fd8
# ╠═fc855167-775a-4688-adbc-93625d72c3cd
# ╠═a5f25a95-d00f-4cc9-b93e-f02efc812e9d
