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
	theme(:juno)

	using LombScargle
end

# ╔═╡ 33b862f3-dc6a-46fe-b73e-a7df7af22e92
using JSON3, SHA

# ╔═╡ 7c9ff7d7-7785-4900-bf29-1bad76947fe1
using LaTeXStrings

# ╔═╡ 30d54c72-876b-4d01-b818-b683ffbc400d
using Roots

# ╔═╡ 234e80df-af67-44ad-8f06-3c3d403dcd25
using KernelDensity

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
pgram = lombscargle(points.day, points.K, samples_per_peak = 100)

# ╔═╡ 55c8d8ef-9d4b-4b9c-8838-b91f1f53f8b0
plot(
	freq(pgram),
	power(pgram),
	xlabel = "частота (1/день)",
	title = "Lomb-Scargle periodogram, сигнал K"
)

# ╔═╡ 5b2930bc-2de0-4388-824a-190d1169cbfe
estimated_period = 2findmaxperiod(pgram)[1]

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

# ╔═╡ 2a542b29-8d3e-49cc-bc8c-2116c6d544f8
initial_params = (;
	mass_quotient = 1.38208,
	observer_angle = 0.981023,
	initial_phase = -1.46867,
	σ_common = [0.023195, 0.0272215],
	offset = [46.6894, 48.8942],
)

# ╔═╡ 30a74a85-c431-469c-bf3d-00190db36c56
channels = [
	ChannelParams(
		measurements_t = points.day,
		measurements_y = points.K,
		darkening_function = claret_darkening,
		darkening_coefs_interpolant = K_coefs_interpolant,
		luminocity_function = phoenixK,
		σ_measured = points.K_err,
		σ_common = FlatPos(0.),
	)
	ChannelParams(
		measurements_t = points.day,
		measurements_y = points.J,
		darkening_function = claret_darkening,
		darkening_coefs_interpolant = J_coefs_interpolant,
		luminocity_function = phoenixJ,
		σ_measured = points.J_err,
		σ_common = FlatPos(0.),
	)
]

# ╔═╡ 97fd2129-d706-480c-a97d-9804027d8b40
model_params = ModelParams(
	channels = channels,
	period = estimated_period,
	β = 0.08,
	temperature_at_bottom = 3500.,
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

# ╔═╡ a4fa7c60-3d71-4091-be53-7e02cdbf59f5
(samples.info.stop_time - samples.info.start_time) / length(samples)

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
	xlabel = "q = m_giant / m_dwarf при априорном q⁻¹ ~ Uniform",
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
	#xlim = (-0.265, -0.23)
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
end;

# ╔═╡ 4cead2e8-2caf-486b-8a5d-990815b88ab9
begin
	local q = 0 : 0.01 : 1.05
	PyPlot.plot(q, i.(q), color = "black", linestyle = "dashed")
	PyPlot.gcf()
end

# ╔═╡ fc855167-775a-4688-adbc-93625d72c3cd
mx = maximize(x -> pdf(k, x[1], x[2]), [1., 60.])

# ╔═╡ a5f25a95-d00f-4cc9-b93e-f02efc812e9d
max_point = Optim.maximizer(mx)

# ╔═╡ 58497ce3-5b59-4473-b3c3-9ddea6de1cb9
begin
	PyPlot.gcf().clear()

	PyPlot.gca().set_xlim(0, 1.8)
	PyPlot.gca().set_ylim(45, 72)
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

# ╔═╡ Cell order:
# ╠═36db1462-6dbf-11ee-38c4-05d52e2c894c
# ╠═b8bda58e-9ed2-4da0-a67a-6d5990e7389d
# ╠═e28e8c98-caa0-41c0-bb15-53c6679dda6d
# ╠═55c8d8ef-9d4b-4b9c-8838-b91f1f53f8b0
# ╠═5b2930bc-2de0-4388-824a-190d1169cbfe
# ╠═2fe448f3-1744-4bbb-83e7-290a9214e7c8
# ╠═d9b9b851-cca0-4a40-a627-7dec9c5da6c1
# ╠═275cb92f-d5d1-4fb9-acb3-5d2317f84a2b
# ╠═232ace15-13a9-4afe-9468-d6e54d796470
# ╠═2a542b29-8d3e-49cc-bc8c-2116c6d544f8
# ╠═00044db4-b168-44be-9d39-87d27b7d330d
# ╠═30a74a85-c431-469c-bf3d-00190db36c56
# ╠═97fd2129-d706-480c-a97d-9804027d8b40
# ╠═c88314a3-cd9e-42b2-acee-4d613b1b36e1
# ╠═eda9134f-b918-42f0-bcfc-e0d601eeeaad
# ╠═a4fa7c60-3d71-4091-be53-7e02cdbf59f5
# ╠═33b862f3-dc6a-46fe-b73e-a7df7af22e92
# ╠═a5070b94-48c2-4405-af78-fddd5784161e
# ╠═174cd8b8-1d1c-4141-a170-6f978f5195e1
# ╠═15ae4b29-3d14-4ad1-811d-de973095f25d
# ╠═45422b39-64d5-4a75-b8c0-8ba0011ba089
# ╠═7c9ff7d7-7785-4900-bf29-1bad76947fe1
# ╠═43973ad5-74c8-4bb7-92c2-372e6fb722dd
# ╠═8ec7dad5-b377-43df-8a7f-ad8c9179d7e2
# ╠═9a192e6e-f8ac-48fc-a7c3-7cbc5ebe7554
# ╠═61b8f08b-4252-4ea7-9b27-37771331de77
# ╠═30d54c72-876b-4d01-b818-b683ffbc400d
# ╠═234e80df-af67-44ad-8f06-3c3d403dcd25
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
# ╠═4cead2e8-2caf-486b-8a5d-990815b88ab9
# ╠═3b39b5a4-6cb1-4b80-9d13-730f6a797fd8
# ╠═fc855167-775a-4688-adbc-93625d72c3cd
# ╠═a5f25a95-d00f-4cc9-b93e-f02efc812e9d
