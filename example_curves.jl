### A Pluto.jl notebook ###
# v0.19.29

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
	Plots.theme(:juno)

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
	mass_quotients = 0.01:0.01:10
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
	local m = interpolated_mesh(1.)
	f = viz(domain(m), color = values(m, 0).g)
end

# ╔═╡ 28452494-35d6-443a-8986-8c41cdc239de
begin
	mesh = interpolated_mesh(1)
	mesh = apply_function(mesh, g -> g^0.25, :g, :T)
	mesh = apply_function(mesh, T -> T^4, :T, :L)

	normals = calc_function_on_faces(mesh, normal)
    areas = calc_function_on_faces(mesh, area)
end

# ╔═╡ 1f6a872e-5411-47e7-a642-97e9180af6c7
md"### Кривые блеска при разных параметрах"

# ╔═╡ 05465b1e-1b8e-4f76-9ae8-791e2de3c050
phases = -π/2 : 0.001 : 3π/2

# ╔═╡ 424af3e5-703b-46a1-be45-e87f78714517
params = (;
	mass_quotient = 0.5,
	observer_angle = π/2,
	temperature_at_bottom = 5000,
	β = 0.25,
	interpolated_mesh,
	luminocity_function = black_body_K,
    darkening_function = claret_darkening,
    darkening_coefficients = (1.3113, -1.2998, 1.0144, -0.3272)
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
		star_magnitude(phases; params..., luminocity_function = black_body_K) .+ 0.02,
		label = "∫planck(λ, T) dλ"
	)
	plot!(
		phases,
		star_magnitude(phases; params..., luminocity_function = T_4) .+ 23.256,
		label = "T^4"
	)
	
end

# ╔═╡ 878919d1-66ad-4ff4-9832-2fbb197ac01f
begin
	plot( xlabel = "phase", ylabel = "m", yflip = true)
	plot!(
		phases,
		star_magnitude(phases; params..., luminocity_function = black_body_K, temperature_at_bottom = 2000),
	)
	plot!(
		phases,
		star_magnitude(phases; params..., luminocity_function = black_body_K, temperature_at_bottom = 3000) .+ 1.22,
	)
	plot!(
		phases,
		star_magnitude(phases; params..., luminocity_function = black_body_K, temperature_at_bottom = 4000) .+ 1.89,
	)
	
end

# ╔═╡ 11872494-739a-4ada-8572-733d1756b6a5
begin
	plot(title = "При разном соотношении масс", xlabel = "phase", ylabel = "m", yflip = true)
	for q ∈ [0.1, 0.2, 0.4, 0.8, 1.6, 3.2]
		plot!(
			phases,
			star_magnitude(phases; params..., mass_quotient = q),
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

# ╔═╡ 8b8880ed-eec4-4ae3-8055-5b899915dff7
begin
	plot(title = "Зависимость светимости от температуры", xlabel = "T", ylabel = "L", legend = :bottomright)
	plot!(temperature_nodes, black_body_K_rectangle.(temperature_nodes), label = "planck(λ, T) * width")
	plot!(temperature_nodes, black_body_K.(temperature_nodes), label = "∫planck(λ, T) dλ")
end

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
	plot!(
		phases,
		star_magnitude(phases; params..., darkening_function = one, darkening_coefficients = ()),
		label = "Без потемнения к краю"
	)
	plot!(
		phases,
		star_magnitude(phases; params...) .- 0.114 ,
		label = "С потемнением T = 3600"
	)
	plot!(
		phases,
		star_magnitude(phases; params...,
			darkening_coefficients = darkening_coefficients2
		) .- 0.19,
		label = "С потемнением T = 2000"
	)
	plot!(legend = :bottomright)
end

# ╔═╡ 6fbfd8b0-19bc-417f-8fb8-9345086685f3
md"### Скорость интегрирования с кешированием нормалей"

# ╔═╡ f457d950-a90f-4726-b132-f5f6ceaee8e3
d = (1/√2, 1/√2, 0.)

# ╔═╡ 9ec17918-f502-435e-ac1a-d7534c9362f3
# Без предварительно вычисленных нормалей
@btime integrate_data_over_triangular_mesh(mesh, :g, d)

# ╔═╡ 8f019dcc-fff8-4a7f-a15c-d5e5d507656f
# С предварительно вычисленными нормалями
@btime integrate_data_over_mesh(mesh, :g, d, normals, areas, c -> 1., ())

# ╔═╡ 16bbab2a-43ff-432a-af10-d83ac1a7aa4b
# С предварительно вычисленными нормалями и учетом потемнения к краю
@btime integrate_data_over_mesh(mesh, :g, d, normals, areas, claret_darkening, darkening_coefficients)

# ╔═╡ c6c9f467-d821-47e2-b2a1-fedfdb685880
@code_warntype integrate_data_over_mesh(mesh, :g, d, normals, areas, claret_darkening, darkening_coefficients)

# ╔═╡ Cell order:
# ╠═0f19eafc-6338-11ee-346c-d781d36c948a
# ╠═243eae5b-8f0a-4246-8f2e-3de35bb3941a
# ╟─643330fc-6808-402c-b204-ba5cc88f378a
# ╠═74d8becd-7529-47dd-bf4b-9b1cfe3f17de
# ╠═87a082a1-4e29-4702-acb0-35afc8e51735
# ╠═c72dec8c-9841-45af-bb0f-c0818102fe4f
# ╠═3a8746b9-7348-4e60-8fa7-ec1ab2fe3841
# ╠═28452494-35d6-443a-8986-8c41cdc239de
# ╟─1f6a872e-5411-47e7-a642-97e9180af6c7
# ╠═05465b1e-1b8e-4f76-9ae8-791e2de3c050
# ╠═424af3e5-703b-46a1-be45-e87f78714517
# ╠═08d08b9d-f632-49f6-8df6-c39d35b7dc27
# ╠═878919d1-66ad-4ff4-9832-2fbb197ac01f
# ╠═11872494-739a-4ada-8572-733d1756b6a5
# ╠═41327e56-5682-438b-bb91-a7acc918d589
# ╠═ae5f7182-0a38-4178-bf15-a7307876826e
# ╠═8b8880ed-eec4-4ae3-8055-5b899915dff7
# ╠═5b03df78-85f1-4f54-b814-9658f4666779
# ╠═d88a8c9a-9afb-4f0e-95d8-3b60f0fcf22e
# ╠═21085b08-dad7-4d62-963e-5355576be10b
# ╠═e1cb3cb6-532d-43f8-beaa-646f8c1264fa
# ╠═4b1e0d61-d58a-4e64-a120-1da21c8ece50
# ╟─6fbfd8b0-19bc-417f-8fb8-9345086685f3
# ╠═f457d950-a90f-4726-b132-f5f6ceaee8e3
# ╠═9ec17918-f502-435e-ac1a-d7534c9362f3
# ╠═8f019dcc-fff8-4a7f-a15c-d5e5d507656f
# ╠═16bbab2a-43ff-432a-af10-d83ac1a7aa4b
# ╠═c6c9f467-d821-47e2-b2a1-fedfdb685880
