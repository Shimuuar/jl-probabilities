### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ 9b761aac-0fac-11ef-1966-633f3bbdcbc1
begin
	import Pkg
	Pkg.activate(".")

	using TuringStars
	
	using FITSIO
	using Trapz

	using Plots
	plotlyjs()
end

# ╔═╡ 28a2baf2-d13e-4ee3-998b-f85b941e45b5
fits_path(Teff) = "atmosphere_tables/phoenix_spec_intensities/lte0$Teff-4.00-0.0.PHOENIX-ACES-AGSS-COND-SPECINT-2011.fits"

# ╔═╡ 8e927295-a8dc-4452-819f-ca096b5fd61f
Teff = 2300

# ╔═╡ e7f68f56-9d36-4255-937b-cb7e127981f0
f = FITS(fits_path(Teff))

# ╔═╡ d30ba23e-3697-4f59-9dc7-48caa875a191
read_header(f[1])

# ╔═╡ e6b57405-43fe-4ab8-944e-437ad207d226
cosines = read(f[2]) # cos θ, последний элемент равен 1

# ╔═╡ c91f8297-b21b-4d8f-93bc-b2608edc571a
intensities = read(f[1])

# ╔═╡ e53a9909-b279-4229-b9aa-f0ff79a91bb9
lambdas = 500 : 500 + size(intensities)[1] - 1

# ╔═╡ 6daf44d8-d9f9-42b2-9f38-7d9e32fb0658
plot(lambdas, intensities[:, end], legend = false)

# ╔═╡ 1783a0f4-5a99-4d8a-a8b0-0af5e39d2393
K_range = findfirst(lambdas .≥ 20_000) : findlast(lambdas .≤ 24_000)

# ╔═╡ 8f78f077-af11-4f14-a26e-c37af1193daf
begin
	local intensitiesK = trapz((lambdas[K_range], :), intensities[K_range, :])
	local normalized_intensities = intensitiesK ./ intensitiesK[end]

	plot(cosines, normalized_intensities, label = "Goettingen")
	plot!(cosines, claret_darkening.(cosines, K_coefs_interpolant(2300)...), label = "Claret")

	plot!(title = "Потемнение к краю при T = $Teff на канале K", xlabel = "cos θ", legend = :bottomright)
end

# ╔═╡ 354902aa-f196-468c-b275-6c36913f2ec4
J_range = findfirst(lambdas .≥ 10_850) : findlast(lambdas .≤ 14_150)

# ╔═╡ b1a88768-566e-482d-8095-e95e60d0fcb7
begin
	Ts = 2300 : 100 : 4900
	Ks = []
	Js = []

	for T ∈ Ts
		file = FITS(fits_path(T))
		spectre = read(file[1])[:, end]

		K = trapz(lambdas[K_range], spectre[K_range])
		push!(Ks, K)
		J = trapz(lambdas[J_range], spectre[J_range])
		push!(Js, J)
	end
end

# ╔═╡ 06507e37-bc8d-43bb-8f6a-39750944311f
plot(
	Ts,
	[Ks, Js],
	label = ["K" "J"],
	title = "L = f(T)",
	legend = :topleft
)

# ╔═╡ bc0b8ef1-d400-4662-91fe-af97a4b98b4e
print(join(string.(Ks), ", "))

# ╔═╡ 17e920f2-4881-494b-ba03-1f208a3963d9
print(join(string.(Js), ", "))

# ╔═╡ Cell order:
# ╠═9b761aac-0fac-11ef-1966-633f3bbdcbc1
# ╠═28a2baf2-d13e-4ee3-998b-f85b941e45b5
# ╠═8e927295-a8dc-4452-819f-ca096b5fd61f
# ╠═e7f68f56-9d36-4255-937b-cb7e127981f0
# ╠═d30ba23e-3697-4f59-9dc7-48caa875a191
# ╠═e6b57405-43fe-4ab8-944e-437ad207d226
# ╠═c91f8297-b21b-4d8f-93bc-b2608edc571a
# ╠═e53a9909-b279-4229-b9aa-f0ff79a91bb9
# ╠═6daf44d8-d9f9-42b2-9f38-7d9e32fb0658
# ╠═8f78f077-af11-4f14-a26e-c37af1193daf
# ╠═1783a0f4-5a99-4d8a-a8b0-0af5e39d2393
# ╠═354902aa-f196-468c-b275-6c36913f2ec4
# ╠═b1a88768-566e-482d-8095-e95e60d0fcb7
# ╠═06507e37-bc8d-43bb-8f6a-39750944311f
# ╠═bc0b8ef1-d400-4662-91fe-af97a4b98b4e
# ╠═17e920f2-4881-494b-ba03-1f208a3963d9
