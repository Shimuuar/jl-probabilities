### A Pluto.jl notebook ###
# v0.19.12

using Markdown
using InteractiveUtils

# ╔═╡ 4506d47e-7c17-46f2-a153-ed3b0b736d11
expr = :(12x+4)

# ╔═╡ be0ec276-16d6-4501-b600-0f8ef2168ace
vars = (;x = 2)

# ╔═╡ d912cc35-c1cc-42f3-b082-2e0242bdf2e9
md"### Вычислитель"

# ╔═╡ aa41b155-bcf5-4aaa-90dc-380023e2ca57
function substitute_vars!(args::Array, vars::NamedTuple)
	for i ∈ 2 : length(args)
		if isa(args[i], Expr)
			substitute_vars!(args[i].args, vars)

		elseif isa(args[i], Symbol)
			try
				args[i] = getfield(vars, args[i])
			catch e
				if !isa(e, KeyError)
					throw(e)
				end
			end
		end
	end
	return args
end

# ╔═╡ 062fee94-160a-4092-ad2d-490582f2110f
function substitute_vars(expr::Expr, vars::NamedTuple)
	expr = deepcopy(expr)
	substitute_vars!(expr.args, vars)
	return expr
end

# ╔═╡ a32df432-35cf-49d3-be88-f65ae7601c7f
substitute_vars(:(12x+4), (; x = 2))

# ╔═╡ e6350425-f8e5-42ac-9aa8-b86ec0dde470
function fold_commuting_consts!(args::Array)
	for i ∈ 2 : length(args)
		if isa(args[i], Expr)
			res = fold_commuting_consts!(args[i].args)
			if res != nothing
				args[i] = res
			end
		end
	end

	where_numbers = isa.(args, Number)
	n_numbers = sum(where_numbers)
	if n_numbers > 1
		func = getfield(Base, args[1])
		res = func(args[where_numbers]...)
		if n_numbers == length(args) - 1
			return res
		else
			vars_and_exprs = args[.!where_numbers]
			resize!(args, length(vars_and_exprs) + 1)
			args[2] = res
			@. args[3:end] = vars_and_exprs[2:end]
			return nothing
		end
	end
end

# ╔═╡ 4c95167f-429b-410e-b777-1925ce3b9b29
function fold_commuting_consts(expr::Expr)
	expr = deepcopy(expr)
	res = fold_commuting_consts!(expr.args)
	if res != nothing
		return res
	else
		return expr
	end
end

# ╔═╡ 5064e19c-ae48-45fd-a1b5-488064e97f01
function evaluate(expr::Expr, vars)
	expr = substitute_vars(expr, vars)
	res = fold_commuting_consts!(expr.args)
	if res != nothing
		return res
	else
		return expr
	end
end

# ╔═╡ 89a6a39f-bf6a-4c4c-9626-b239dc70bbb5
fold_commuting_consts(:(x + 2 + 3*18*y + 4))

# ╔═╡ f87cec96-b7ee-4dfa-a6d5-2ce508857bce
evaluate(:(12x+4), (; x = 2))

# ╔═╡ 75296fe5-db59-44b2-885a-9c2894757f8b
evaluate(:(x + 2 + 3*18*y + 4), (; x = 0.5, y = -1))

# ╔═╡ d0b16392-8079-4915-bdf9-00fcbd6c0f2f
0.5 + 2 + 3*18*(-1) + 4

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.2"
manifest_format = "2.0"
project_hash = "da39a3ee5e6b4b0d3255bfef95601890afd80709"

[deps]
"""

# ╔═╡ Cell order:
# ╠═4506d47e-7c17-46f2-a153-ed3b0b736d11
# ╠═be0ec276-16d6-4501-b600-0f8ef2168ace
# ╟─d912cc35-c1cc-42f3-b082-2e0242bdf2e9
# ╠═aa41b155-bcf5-4aaa-90dc-380023e2ca57
# ╠═062fee94-160a-4092-ad2d-490582f2110f
# ╠═a32df432-35cf-49d3-be88-f65ae7601c7f
# ╠═e6350425-f8e5-42ac-9aa8-b86ec0dde470
# ╠═4c95167f-429b-410e-b777-1925ce3b9b29
# ╠═5064e19c-ae48-45fd-a1b5-488064e97f01
# ╠═89a6a39f-bf6a-4c4c-9626-b239dc70bbb5
# ╠═f87cec96-b7ee-4dfa-a6d5-2ce508857bce
# ╠═75296fe5-db59-44b2-885a-9c2894757f8b
# ╠═d0b16392-8079-4915-bdf9-00fcbd6c0f2f
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
