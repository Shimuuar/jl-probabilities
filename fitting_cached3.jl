### A Pluto.jl notebook ###
# v0.19.36

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
	theme(:default)

	using LombScargle
end

# ╔═╡ 33b862f3-dc6a-46fe-b73e-a7df7af22e92
using JSON3, SHA

# ╔═╡ 234e80df-af67-44ad-8f06-3c3d403dcd25
using KernelDensity

# ╔═╡ 360119a3-9cba-4e40-ad6c-88a0be627b1e
using Roots

# ╔═╡ 3b39b5a4-6cb1-4b80-9d13-730f6a797fd8
using Optim

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

# ╔═╡ 7e8b804f-b511-4359-8c44-286743a2ff9b
estimated_period

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
				darkening_coefficients = channel.darkening_coefficients
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
		σ_measured = points.J_err,
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
	n_samples = 4096,
	init_params = initial_params,
	sampler = NUTS()
)

# ╔═╡ eda9134f-b918-42f0-bcfc-e0d601eeeaad
samples = cached_sample(chain_params)

# ╔═╡ a5070b94-48c2-4405-af78-fddd5784161e
chain_params |> JSON3.write |> sha1 |> bytes2hex

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
	xlim = (-0.265, -0.23)
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
	samples[:mass_quotient],
	samples[:observer_angle],
	markersize = 0.2,
	xlabel = "Масса карлика / масса гиганта",
	ylabel = "Наклонение (рад.)",
	legend = false
)

# ╔═╡ fcff208a-6eea-475b-8cee-908457bbc1d3
k = kde((
	1 ./ samples[:mass_quotient].data[:, 1],
	samples[:observer_angle].data[:, 1] .* 180 ./ π
))

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
labels_ = [0.95, 0.68, 0]#, 0.38]

# ╔═╡ 230396d1-c1ce-4b71-a533-ceb745392393
text_lavels = @. string(Int(labels_ * 100)) * "%"

# ╔═╡ f5bccc4e-aefe-4bc7-a6c2-da3de047b5d4
thresholds = get_threshold.(Ref(k), labels_)

# ╔═╡ 5d7caa39-2328-4e12-805b-3438ce74faa8
import PyPlot

# ╔═╡ 072d2ba5-4bf8-4d9d-8165-38bcbfe8cd82
PyPlot.svg(true)

# ╔═╡ 4693a488-e017-46a8-b077-cffd3b20303e
function i(q; m = 1.44)
	rad2deg(asin(cbrt(0.322408 / m * (1 + q)^2)))
end

# ╔═╡ 02435262-6c25-427f-9de6-b5251067eccb
samples_filtered = filter(DataFrame(samples)) do sample
	try
		i_min = i(1 / sample[:mass_quotient])
		return i_min < rad2deg(sample[:observer_angle])
	catch DomainError
		return false
	end
end

# ╔═╡ d6d41a41-6528-49b2-bf52-3f80542f9dc0
begin
	local q = 0 : 0.01 : 1.05
	plot(q, i.(q))
	plot!(q, i.(q, m = 1.41))
end

# ╔═╡ 439fba04-2fd8-47ef-b34c-d0b63d9508ce
i(1.113)

# ╔═╡ 4cead2e8-2caf-486b-8a5d-990815b88ab9
begin
	local q = 0 : 0.01 : 1.05
	PyPlot.plot(q, i.(q), color = "black", linestyle = "dashed")
	PyPlot.gcf()
end

# ╔═╡ 7eb4f613-beed-47a3-995a-c16b7219e597
Dict(zip(
		thresholds,
		text_lavels
	))

# ╔═╡ 76cc35b4-d81c-4ab6-ac3e-c0f4ada37de8


# ╔═╡ fc855167-775a-4688-adbc-93625d72c3cd
mx = maximize(x -> pdf(k, x[1], x[2]), [1., 60.])

# ╔═╡ a5f25a95-d00f-4cc9-b93e-f02efc812e9d
max_point = Optim.maximizer(mx)

# ╔═╡ 58497ce3-5b59-4473-b3c3-9ddea6de1cb9
begin
	PyPlot.gcf().clear()

	PyPlot.gca().set_xlim(0, 1.8)
	PyPlot.gca().set_ylim(42, 70)
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

# ╔═╡ ad9d27f1-3fce-4cb1-bcb4-487d3f91c344
k_filtered = kde((
	1 ./ samples_filtered.mass_quotient,
	samples_filtered.observer_angle .* 180 ./ π
))

# ╔═╡ 9369dee7-b21f-4afc-929b-7633fe80f3af
mx_filtered = maximize(x -> pdf(k_filtered, x[1], x[2]), [0.5, 50.])

# ╔═╡ 65375d07-0d16-454f-97d7-822aa55e9c84
max_point_filtered = Optim.maximizer(mx_filtered)

# ╔═╡ b34c8f5d-c125-4ebf-aec0-9b37d5432dc8
thresholds_filtered = get_threshold.(Ref(k_filtered), labels_)

# ╔═╡ af369c17-81b1-47bc-ae3d-af7212a34afb
begin
	PyPlot.gcf().clear()

	# PyPlot.gca().set_xlim(0, 1.8)
	# PyPlot.gca().set_ylim(42, 70)
	PyPlot.gca().set_xlabel("Масса гиганта / масса карлика")
	PyPlot.gca().set_ylabel("Наклонение (°)")

	local cs = PyPlot.contour(
		k_filtered.x, k_filtered.y, k_filtered.density',
		levels = thresholds_filtered
	)
	PyPlot.clabel(cs, fmt = Dict(zip(
		thresholds_filtered,
		text_lavels
	)))
	PyPlot.scatter([max_point_filtered[1]], [max_point_filtered[2]])

	local q = 0 : 0.01 : 0.8
	PyPlot.plot(q, i.(q))
	PyPlot.gcf()
end

# ╔═╡ b5cf0296-5308-42f0-b9b4-0f2d7ba0c0c0
begin
	scatter(
		1 ./ samples_filtered.mass_quotient,
		rad2deg.(samples_filtered.observer_angle),
	)
	local q = 0.1 : 0.01 : 0.7
	plot!(q, i.(q), legend = false)
end

# ╔═╡ 60edb60f-8951-4d92-8395-c99d808b9e29
density(1 ./ samples_filtered.mass_quotient, xlabel = "Масса гиганта / масса карлика")

# ╔═╡ 2c22988c-b0c4-4003-a1d9-5f902d31c980
dwarf_mass(q, i) = 0.322408 * (1 + q)^2 / sin(i)^3

# ╔═╡ 2a7073d9-afb1-4d08-a76e-fb88248e4a26
begin
	dwarf_masses = dwarf_mass.(1 ./ samples_filtered.mass_quotient, samples_filtered.observer_angle)

	giant_masses = dwarf_masses ./ samples_filtered.mass_quotient

	scatter(dwarf_masses, giant_masses, markersize = 1)

	#plot!(xlim = (0., 1.45), ylim = (0, 1))
end

# ╔═╡ Cell order:
# ╠═36db1462-6dbf-11ee-38c4-05d52e2c894c
# ╠═b8bda58e-9ed2-4da0-a67a-6d5990e7389d
# ╠═e28e8c98-caa0-41c0-bb15-53c6679dda6d
# ╠═55c8d8ef-9d4b-4b9c-8838-b91f1f53f8b0
# ╠═5b2930bc-2de0-4388-824a-190d1169cbfe
# ╠═7e8b804f-b511-4359-8c44-286743a2ff9b
# ╠═2fe448f3-1744-4bbb-83e7-290a9214e7c8
# ╠═d9b9b851-cca0-4a40-a627-7dec9c5da6c1
# ╠═275cb92f-d5d1-4fb9-acb3-5d2317f84a2b
# ╠═232ace15-13a9-4afe-9468-d6e54d796470
# ╠═960ab30d-a1fa-4803-a4d4-d0860286ba87
# ╠═00044db4-b168-44be-9d39-87d27b7d330d
# ╠═30a74a85-c431-469c-bf3d-00190db36c56
# ╠═97fd2129-d706-480c-a97d-9804027d8b40
# ╠═c88314a3-cd9e-42b2-acee-4d613b1b36e1
# ╠═eda9134f-b918-42f0-bcfc-e0d601eeeaad
# ╠═33b862f3-dc6a-46fe-b73e-a7df7af22e92
# ╠═a5070b94-48c2-4405-af78-fddd5784161e
# ╠═174cd8b8-1d1c-4141-a170-6f978f5195e1
# ╠═15ae4b29-3d14-4ad1-811d-de973095f25d
# ╠═45422b39-64d5-4a75-b8c0-8ba0011ba089
# ╠═43973ad5-74c8-4bb7-92c2-372e6fb722dd
# ╠═8ec7dad5-b377-43df-8a7f-ad8c9179d7e2
# ╠═9a192e6e-f8ac-48fc-a7c3-7cbc5ebe7554
# ╠═61b8f08b-4252-4ea7-9b27-37771331de77
# ╠═234e80df-af67-44ad-8f06-3c3d403dcd25
# ╠═360119a3-9cba-4e40-ad6c-88a0be627b1e
# ╠═02435262-6c25-427f-9de6-b5251067eccb
# ╠═fcff208a-6eea-475b-8cee-908457bbc1d3
# ╠═46a12905-4167-4f7e-92a6-6a633903f7e9
# ╠═dedc4007-cf2e-49ac-860f-b6b27cd2be05
# ╠═230396d1-c1ce-4b71-a533-ceb745392393
# ╠═f5bccc4e-aefe-4bc7-a6c2-da3de047b5d4
# ╠═5d7caa39-2328-4e12-805b-3438ce74faa8
# ╠═072d2ba5-4bf8-4d9d-8165-38bcbfe8cd82
# ╠═58497ce3-5b59-4473-b3c3-9ddea6de1cb9
# ╠═4693a488-e017-46a8-b077-cffd3b20303e
# ╠═d6d41a41-6528-49b2-bf52-3f80542f9dc0
# ╠═439fba04-2fd8-47ef-b34c-d0b63d9508ce
# ╠═4cead2e8-2caf-486b-8a5d-990815b88ab9
# ╠═7eb4f613-beed-47a3-995a-c16b7219e597
# ╠═76cc35b4-d81c-4ab6-ac3e-c0f4ada37de8
# ╠═3b39b5a4-6cb1-4b80-9d13-730f6a797fd8
# ╠═fc855167-775a-4688-adbc-93625d72c3cd
# ╠═9369dee7-b21f-4afc-929b-7633fe80f3af
# ╠═a5f25a95-d00f-4cc9-b93e-f02efc812e9d
# ╠═65375d07-0d16-454f-97d7-822aa55e9c84
# ╠═ad9d27f1-3fce-4cb1-bcb4-487d3f91c344
# ╠═b34c8f5d-c125-4ebf-aec0-9b37d5432dc8
# ╠═af369c17-81b1-47bc-ae3d-af7212a34afb
# ╠═b5cf0296-5308-42f0-b9b4-0f2d7ba0c0c0
# ╠═60edb60f-8951-4d92-8395-c99d808b9e29
# ╠═2c22988c-b0c4-4003-a1d9-5f902d31c980
# ╠═2a7073d9-afb1-4d08-a76e-fb88248e4a26
