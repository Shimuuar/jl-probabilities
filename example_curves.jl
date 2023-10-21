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
	local m = interpolated_mesh(1)
	f = viz(domain(m), color = values(m, 0).g)
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
	luminocity_function = T_4
)

# ╔═╡ 08d08b9d-f632-49f6-8df6-c39d35b7dc27
begin
	plot(title = "При разной функции L = f(T)", xlabel = "phase", ylabel = "m", yflip = true)
	plot!(
		phases,
		star_magnitude(phases; params..., luminocity_function = black_body_K_rectangle),
		label = "planck(T) * width"
	)
	plot!(
		phases,
		star_magnitude(phases; params..., luminocity_function = T_4) .+ 23.266,
		label = "T^4"
	)
	
end

# ╔═╡ 11872494-739a-4ada-8572-733d1756b6a5
begin
	plot(title = "При разном соотношении масс; L = T^4", xlabel = "phase", ylabel = "m", yflip = true)
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
	plot(title = "При разном наклонении; L = T^4", xlabel = "phase", ylabel = "m", yflip = true)
	for i ∈ [0., π/6, π/4, π/3, π/2]
		plot!(
			phases,
			star_magnitude(phases; params..., mass_quotient = 1, observer_angle = i),
			label = "i = $(round(i, digits = 3))"
		)
	end
	plot!()
end

# ╔═╡ f17ee0d2-95dc-4662-8ce9-904c5ef56b1f


# ╔═╡ Cell order:
# ╠═0f19eafc-6338-11ee-346c-d781d36c948a
# ╠═243eae5b-8f0a-4246-8f2e-3de35bb3941a
# ╟─643330fc-6808-402c-b204-ba5cc88f378a
# ╠═74d8becd-7529-47dd-bf4b-9b1cfe3f17de
# ╠═87a082a1-4e29-4702-acb0-35afc8e51735
# ╠═c72dec8c-9841-45af-bb0f-c0818102fe4f
# ╠═3a8746b9-7348-4e60-8fa7-ec1ab2fe3841
# ╟─1f6a872e-5411-47e7-a642-97e9180af6c7
# ╠═05465b1e-1b8e-4f76-9ae8-791e2de3c050
# ╠═424af3e5-703b-46a1-be45-e87f78714517
# ╠═08d08b9d-f632-49f6-8df6-c39d35b7dc27
# ╠═11872494-739a-4ada-8572-733d1756b6a5
# ╠═41327e56-5682-438b-bb91-a7acc918d589
# ╠═f17ee0d2-95dc-4662-8ce9-904c5ef56b1f
