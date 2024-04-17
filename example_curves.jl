### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ 0f19eafc-6338-11ee-346c-d781d36c948a
begin
	import Pkg
	Pkg.activate(".")

	using Revise
	using TuringStars

	using Plots
	using StatsPlots
	plotlyjs()
	Plots.theme(:default)

	using Meshes
	using GeoTables

	import Makie
	import WGLMakie

	using BenchmarkTools
end

# ╔═╡ 243eae5b-8f0a-4246-8f2e-3de35bb3941a
using LinearAlgebra

# ╔═╡ 643330fc-6808-402c-b204-ba5cc88f378a
md"### Визуализация направлений"

# ╔═╡ 74d8becd-7529-47dd-bf4b-9b1cfe3f17de
begin
	mass_quotients = 0.001:0.01:10
	interpolated_mesh = InterpolatedRocheMesh(64, mass_quotients)
end

# ╔═╡ 87a082a1-4e29-4702-acb0-35afc8e51735
begin
	local observer_angle = pi/4
	local phases = 0 : 0.2 : 2*π
	directions = [(
		sin(observer_angle) * cos(phase),
		sin(observer_angle) * sin(phase),
		cos(observer_angle)
	) for phase ∈ phases]
end

# ╔═╡ c72dec8c-9841-45af-bb0f-c0818102fe4f
begin
	local f = viz(domain(interpolated_mesh(0.5)), showfacets = true)

	for d in directions
		Makie.lines!([d, (0, 0, 0)])
	end

	f
end

# ╔═╡ 3a8746b9-7348-4e60-8fa7-ec1ab2fe3841
begin
	m = interpolated_mesh(1.)
	m = avg_over_faces(m, :g)
	m = apply_function(m, g -> g^0.08, :g, :T)
	f = viz(domain(m), color = values(m, 2).T, showfacets = true)
	# cbar(f, values(m, 0).T)
end

# ╔═╡ 692e28c0-ca1a-4186-8009-342a2e4efdb0
begin
	local temperatures = values(m, 2).T .* 3500.
	plot(
		temperatures,
		zcolor = temperatures,
		colormap = :viridis,
		format = :svg
	)
end

# ╔═╡ 1f6a872e-5411-47e7-a642-97e9180af6c7
md"### Кривые блеска при разных параметрах"

# ╔═╡ 05465b1e-1b8e-4f76-9ae8-791e2de3c050
phases = -π/2 : 0.001 : 3π/2

# ╔═╡ 424af3e5-703b-46a1-be45-e87f78714517
params = (;
	mass_quotient = 0.5,
	observer_angle = π/3,
	temperature_at_bottom = 3500,
	β = 0.08,
	interpolated_mesh,
	luminocity_function = black_body_K,
    darkening_function = claret_darkening,
	# darkening_coefs_interpolant = K_coefs_interpolant,
	darkening_coefs_interpolant = T -> (1.3113, -1.2998, 1.0144, -0.3272),
)

# ╔═╡ 08d08b9d-f632-49f6-8df6-c39d35b7dc27
begin
	plot(title = "При разной функции L = f(T)", xlabel = "phase", ylabel = "m", yflip = true)
	plot!(
		phases,
		star_magnitude(phases; params..., luminocity_function = black_body_K_rectangle),
		label = "planck(λ, T) * width"
	)
	plot!(
		phases,
		star_magnitude(phases; params..., luminocity_function = black_body_K) .+ 0.015,
		label = "∫planck(λ, T) dλ"
	)
	plot!(
		phases,
		star_magnitude(phases; params..., luminocity_function = T_4) .+ 22.4,
		label = "T^4"
	)
	
end

# ╔═╡ 878919d1-66ad-4ff4-9832-2fbb197ac01f
begin
	plot(xlabel = "phase", ylabel = "m", yflip = true)
	magnitudes = star_magnitude(phases; params..., temperature_at_bottom = 2000)
	plot!(
		phases,
		magnitudes .- magnitudes[1],
		label = "T_bottom = 2000"
	)

	magnitudes = star_magnitude(phases; params..., temperature_at_bottom = 3000)
	plot!(
		phases,
		magnitudes .- magnitudes[1],
		label = "T_bottom = 3000"
	)

	magnitudes = star_magnitude(phases; params..., temperature_at_bottom = 4000)
	plot!(
		phases,
		magnitudes .- magnitudes[1],
		label = "T_bottom = 4000"
	)
	
end

# ╔═╡ 11872494-739a-4ada-8572-733d1756b6a5
begin
	plot(title = "При разном соотношении масс", xlabel = "phase", ylabel = "m", yflip = true)
	for q ∈ [0.1, 0.2, 0.4, 0.8, 1.6, 3.2]
		plot!(
			phases,
			star_magnitude(phases; params..., mass_quotient = q, observer_angle = pi/3),
			label = "q = $q"
		)
	end
	plot!()
end

# ╔═╡ 41327e56-5682-438b-bb91-a7acc918d589
begin
	plot(title = "При разном наклонении", xlabel = "phase", ylabel = "m", yflip = true)
	for i ∈ [0., π/6, π/4, π/3, π/2]
		plot!(
			phases,
			star_magnitude(phases; params..., mass_quotient = 1, observer_angle = i),
			label = "i = $(round(i, digits = 3))"
		)
	end
	plot!()
end

# ╔═╡ ae5f7182-0a38-4178-bf15-a7307876826e
temperature_nodes = 0 : 100 : 50_000

# ╔═╡ 5b03df78-85f1-4f54-b814-9658f4666779
darkening_coefficients = (1.3113, -1.2998, 1.0144, -0.3272)

# ╔═╡ d88a8c9a-9afb-4f0e-95d8-3b60f0fcf22e
darkening_coefficients2 = (1.2143, -0.3175, 0.0665, -0.0100)
	# (1.2290, -0.2192, -0.1071, 0.0679)

# ╔═╡ 21085b08-dad7-4d62-963e-5355576be10b
darkening_coefficients3 = (1.6389, -2.4779, 2.1452, -0.7306)

# ╔═╡ e1cb3cb6-532d-43f8-beaa-646f8c1264fa
begin
	local θ = 0 : 0.01 : π/2
	local cosines = cos.(θ)
	plot(xlabel = "θ", title = "Коэффициент потемнения к краю", legend = :bottomleft)
	plot!(θ, claret_darkening.(cosines, darkening_coefficients...), label = "T = 3600, log g = 4.0")
	plot!(θ, claret_darkening.(cosines, darkening_coefficients2...), label = "T = 2000, log g = 4.0")
	plot!(θ, claret_darkening.(cosines, darkening_coefficients3...), label = "T = 3600, log g = 5.0")
end

# ╔═╡ 4b1e0d61-d58a-4e64-a120-1da21c8ece50
begin
	plot(title = "С потемнением к краю и без", xlabel = "phase", ylabel = "m", yflip = true)
	mgs = star_magnitude(phases; params..., darkening_function = one, darkening_coefs_interpolant = T -> ())
	plot!(
		phases,
		mgs .- mgs[1],
		label = "Без потемнения к краю"
	)

	mgs = star_magnitude(phases; params..., darkening_coefs_interpolant = T -> darkening_coefficients)
	plot!(
		phases,
		mgs .- mgs[1],
		label = "С потемнением T = 3600"
	)

	mgs = star_magnitude(phases; params...,
			darkening_coefs_interpolant = T -> darkening_coefficients2
		)
	plot!(
		phases,
		mgs .- mgs[1],
		label = "С потемнением T = 2000"
	)

	 mgs = star_magnitude(phases; params...)
	plot!(
		phases,
		mgs .- mgs[1],
		label = "С интерполированным потемнением"
	)

	mgs = star_magnitude(phases; params...,
			darkening_coefs_interpolant = T -> K_coefs_interpolant(2500.)
		)
	plot!(
		phases,
		mgs .- mgs[1],
		label = "С потемнением T = 2500"
	)
	plot!(legend = :bottomright)
end

# ╔═╡ a5d7d64d-e2a4-49bc-9541-543acf527403
begin
	plot(title = "С потемнением к краю и без", xlabel = "phase", ylabel = "m", yflip = true)
	local mgs = star_magnitude(phases; params..., darkening_function = one, darkening_coefs_interpolant = T -> ())
	plot!(
		phases,
		mgs .- mgs[1],
		label = "Без потемнения к краю"
	)

	for T in 2000 : 200 : 3600
		darkening_coefficients = K_coefs_interpolant(T)
		mgs = star_magnitude(phases; params..., darkening_coefs_interpolant = T -> darkening_coefficients)
	
		plot!(
			phases,
			mgs .- mgs[1],
			label = "T = $T"
		)

	end

	mgs = star_magnitude(phases; params...)
	plot!(
		phases,
		mgs .- mgs[1],
		label = "С интерполированным потемнением"
	)
	plot!(legend = false)
end

# ╔═╡ 6fbfd8b0-19bc-417f-8fb8-9345086685f3
md"### Скорость интегрирования с кешированием нормалей"

# ╔═╡ 2c716424-cf65-4141-b9ed-e282698ec438
begin
	normals = calc_function_on_faces(m, normalized_normal)
    areas = calc_function_on_faces(m, area)
	mesh = apply_function(m, K_coefs_interpolant, :T, :darkening_coefs)
end

# ╔═╡ f457d950-a90f-4726-b132-f5f6ceaee8e3
d = (1/√2, 1/√2, 0.)

# ╔═╡ 16bbab2a-43ff-432a-af10-d83ac1a7aa4b
@btime integrate_data_over_mesh(mesh, :g, d, normals, areas, claret_darkening)

# ╔═╡ 8f373d53-15f4-4c6a-91fe-cc12c9726285
@code_warntype integrate_data_over_mesh(mesh, :g, d, normals, areas, claret_darkening, K_coefs_interpolant)

# ╔═╡ e84edaed-0c3b-4f2c-9f6e-e7d1880d0d31
@btime star_magnitude(phases; params...)

# ╔═╡ Cell order:
# ╠═0f19eafc-6338-11ee-346c-d781d36c948a
# ╠═243eae5b-8f0a-4246-8f2e-3de35bb3941a
# ╟─643330fc-6808-402c-b204-ba5cc88f378a
# ╠═74d8becd-7529-47dd-bf4b-9b1cfe3f17de
# ╠═87a082a1-4e29-4702-acb0-35afc8e51735
# ╠═c72dec8c-9841-45af-bb0f-c0818102fe4f
# ╠═3a8746b9-7348-4e60-8fa7-ec1ab2fe3841
# ╠═692e28c0-ca1a-4186-8009-342a2e4efdb0
# ╟─1f6a872e-5411-47e7-a642-97e9180af6c7
# ╠═05465b1e-1b8e-4f76-9ae8-791e2de3c050
# ╠═424af3e5-703b-46a1-be45-e87f78714517
# ╠═08d08b9d-f632-49f6-8df6-c39d35b7dc27
# ╠═878919d1-66ad-4ff4-9832-2fbb197ac01f
# ╠═11872494-739a-4ada-8572-733d1756b6a5
# ╠═41327e56-5682-438b-bb91-a7acc918d589
# ╠═ae5f7182-0a38-4178-bf15-a7307876826e
# ╠═5b03df78-85f1-4f54-b814-9658f4666779
# ╠═d88a8c9a-9afb-4f0e-95d8-3b60f0fcf22e
# ╠═21085b08-dad7-4d62-963e-5355576be10b
# ╠═e1cb3cb6-532d-43f8-beaa-646f8c1264fa
# ╠═4b1e0d61-d58a-4e64-a120-1da21c8ece50
# ╠═a5d7d64d-e2a4-49bc-9541-543acf527403
# ╟─6fbfd8b0-19bc-417f-8fb8-9345086685f3
# ╠═2c716424-cf65-4141-b9ed-e282698ec438
# ╠═f457d950-a90f-4726-b132-f5f6ceaee8e3
# ╠═16bbab2a-43ff-432a-af10-d83ac1a7aa4b
# ╠═8f373d53-15f4-4c6a-91fe-cc12c9726285
# ╠═e84edaed-0c3b-4f2c-9f6e-e7d1880d0d31
