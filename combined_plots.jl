### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ 4a348134-1a98-11ef-290d-933d7a2a61ac
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
	using LaTeXStrings

	using KernelDensity
	import PyPlot
	using Roots
	using Optim
end

# ╔═╡ b1d9141a-c674-4087-aaa6-d372fe0c8270
begin
	points = readdlm("stars/T_CrB_JK.dat")[2:end, :]
	points = map(x -> isa(x, Number) ? x : missing, points)
	points = DataFrame(points, [:day, :J, :J_err, :K, :K_err])
	points = dropmissing(points)
	points.day .-= points.day[1]
	points
end

# ╔═╡ b8772986-d938-4bc4-b71a-4d0a393e8032
period = 227.5687

# ╔═╡ 92329f85-a715-4453-b4da-910b7c6786bc
mesh = InterpolatedRocheMesh(MeshParams());

# ╔═╡ ec0ba187-904b-4253-b7d9-82815cd4af58
model_0 = ModelParams(; period);

# ╔═╡ 1b51f03d-fdef-4a1d-aa6a-2d0a352db9f9
model_q_uniform = ModelParams(; period, model_function = q_uniform_model);

# ╔═╡ 1baeb4e6-184b-4ae6-8fb1-f381651c2f39
model_q_inverted = ModelParams(; period, model_function = q_inverted_model);

# ╔═╡ edf1a809-72e9-4b81-9e73-2abf3658210d
channels = [
	ChannelParams(
		measurements_t = points.day,
		measurements_y = points.K,
		darkening_function = claret_darkening,
		darkening_coefs_interpolant = K_coefs_interpolant,
		luminocity_function = phoenixK,
		σ_measured = points.K_err,
	)
	ChannelParams(
		measurements_t = points.day,
		measurements_y = points.J,
		darkening_function = claret_darkening,
		darkening_coefs_interpolant = J_coefs_interpolant,
		luminocity_function = phoenixJ,
		σ_measured = points.J_err,
	)
];

# ╔═╡ 5fcb704e-eb42-42ea-825b-2f3a572f5c75
samples_0 = load_cache("f8053865657ad849f1830419bdc305135f7a82e1");

# ╔═╡ d1023bb5-0b12-470c-85ca-43a728751671
samples_q_uniform = load_cache("247cfd8ebe8ef323a936a48ed6cdf35d465242fd");

# ╔═╡ c644dbd3-3781-40ba-bd25-93a44a7d3771
samples_q_inv = load_cache("0b8e82061b58a6a038d7039c10097c9951f26d23");

# ╔═╡ 8ca30dfd-18df-4048-bd60-0ceffd802648
begin
	gr()
	mass_ratio_plot = density(title = L"q = m_{\mathrm{giant}} / m_{\mathrm{dwarf}}")
	density!(samples_q_inv[:mass_quotient_inv], label = L"q\ \ \ \ \sim \mathrm{Uniform}")
	density!(samples_q_uniform[:mass_quotient_inv], label = L"q^{-1} \sim \mathrm{Uniform}")
	density!(samples_0[:mass_quotient_inv], label = L"m_{\mathrm{giant}}, m_{\mathrm{dwarf}} \sim \mathrm{Uniform}")
end

# ╔═╡ Cell order:
# ╠═4a348134-1a98-11ef-290d-933d7a2a61ac
# ╠═b1d9141a-c674-4087-aaa6-d372fe0c8270
# ╠═b8772986-d938-4bc4-b71a-4d0a393e8032
# ╠═92329f85-a715-4453-b4da-910b7c6786bc
# ╠═ec0ba187-904b-4253-b7d9-82815cd4af58
# ╠═1b51f03d-fdef-4a1d-aa6a-2d0a352db9f9
# ╠═1baeb4e6-184b-4ae6-8fb1-f381651c2f39
# ╠═edf1a809-72e9-4b81-9e73-2abf3658210d
# ╠═5fcb704e-eb42-42ea-825b-2f3a572f5c75
# ╠═d1023bb5-0b12-470c-85ca-43a728751671
# ╠═c644dbd3-3781-40ba-bd25-93a44a7d3771
# ╠═8ca30dfd-18df-4048-bd60-0ceffd802648
